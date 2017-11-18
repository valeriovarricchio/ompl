/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author :  Valerio Varricchio <valerio@mit.edu> */

#ifndef OMPL_BASE_SPACES_REEDS_SHEPP_MANIFOLD_
#define OMPL_BASE_SPACES_REEDS_SHEPP_MANIFOLD_

#include "ompl/base/spaces/ReedsSheppStateSpace.h"
#include <boost/math/constants/constants.hpp>

namespace ompl
{
    namespace base
    {
        /* _T must be a StateSpace type with a metric defined */
        // TODO can we concept-check at compile time?
        template <class _T>
        class SubRiemannianManifold : public _T // TODO move out into SubRiemannianManifold.h
        {

          public:
            typedef typename _T::StateType StateType;
            typedef std::vector<double> TangentVector;
            typedef std::vector<PrivilegedCoordinate*> coordinates;

            SubRiemannianManifold(const _T& _SS): SS_(_SS){}

            virtual ~SubRiemannianManifold(){} // TODO check this out

            class Box{
                Box(const _T& center, double size);
                virtual bool intersectsHyperplane(const StateType& conf,
                                      const std::vector<double>& normal) = 0;
            };

          protected:
            _T& SS;
            struct PrivilegedCoordinate{
                virtual std::string getName(){ return {"unnamed"}; }
                virtual unsigned int getWeight() = 0;
                virtual TangentVector getTangent(const StateType& center) const = 0;
                /* see if you need the following ... */
                /* virtual double inBoxCoefficient(double t) = 0;
                   virtual double outBoxCoefficient(double t) = 0;
                   virtual double compute(const StateType& center,
                                          const StateType& target) = 0;*/
            };
        };

        class ReedsSheppManifold : public SubRiemannianManifold<ReedsSheppStateSpace>
        {

        private:
            const unsigned int dim_;
            double R_;
            LateralCoordinate l_;
            LongitudinalCoordinate f_;
            HeadingCoordinate h_;

            typedef SubRiemannianManifold<ReedsSheppStateSpace> Base ;
        public:
            SubRiemannianManifold(const ReedsSheppStateSpace& _RS):
                SS(_RS),
                dim_(SS.getDimension()), // could actually write dim_(3) directly
                R_(SS.getMinTurningRadius())
            {
                coordinates.push_back(&f_);
                coordinates.push_back(&h_);
                coordinates.push_back(&l_);
            }

            class ReedsSheppBox : public Base::Box {
                ReedsSheppBox(const _T& center, double size);
                virtual bool intersectsHyperplane(const StateType& center,
                                      const std::vector<double>& normal) = 0;
            };

            struct LongitudinalCoordinate : public Base::PrivilegedCoordinate{
                unsigned int getWeight()
                {
                    return 1;
                }

                std::string getName()
                {
                    return {"longitudinal"};
                }

                TangentVector getTangent(const StateType& center)
                {
                    return {cos(center.getYaw()), sin(center.getYaw()), 0.0};
                }
            };

            struct HeadingCoordinate : public Base::PrivilegedCoordinate{
                unsigned int getWeight()
                {
                    return 1;
                }
                virtual std::string getName(){
                    return {"heading"};
                }

                TangentVector getTangent(const StateType& center)
                {
                    return {0.0, 0.0, 1.0};
                }
            };

            struct LateralCoordinate : public Base::PrivilegedCoordinate{
                unsigned int getWeight()
                {
                    return 2;
                }

                std::string getName()
                {
                    return {"lateral"};
                }

                TangentVector getTangent(const StateType& center)
                {
                    return {-sin(center.getYaw()), cos(center.getYaw()), 0.0};
                }
            };
        };
    } // end of base namespace
} // end of ompl namespace

#endif
