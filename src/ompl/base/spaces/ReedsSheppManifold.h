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

#include <ompl/util/VecUtils.h>

namespace ompl
{
    namespace base
    {
        /* _T must be a StateSpace type with a metric defined */
        // TODO can we concept-check at compile time?
        template <class _T>
        class SubRiemannianManifold : public _T // TODO move out into SubRiemannianManifold.h
        {

          using _T::_T; // inherit base ctors

          protected:
            struct PrivilegedCoordinate;

          public:
            typedef typename ompl::base::State* StatePtr;
            typedef std::vector<double> TangentVector;
            std::vector<PrivilegedCoordinate*> coordinates;

            virtual ~SubRiemannianManifold(){} // TODO check this out

            virtual bool inPositiveHalfspace(const StatePtr state,
                                             const StatePtr pivot,
                                             const TangentVector normal) const = 0;

            virtual TangentVector getSplittingNormal(const StatePtr state,
                                                     uint depth){
                checkSetup();
                return coordinates[lie_split_idx[depth%W_]]->getTangent(state);
            }

            class OuterBox{
              protected:
                double size;
                StatePtr center;
              public:
                OuterBox(const StatePtr center_, double size_):
                    size(size_), center(center_) // TODO note only copying pointer, relying on its validity!
                {}

                virtual bool intersectsHyperplane(const StatePtr state,
                                      const std::vector<double>& normal) const = 0;

                virtual ~OuterBox() = default;
            };

            typedef std::shared_ptr<OuterBox> OuterBoxPtr;

            virtual OuterBoxPtr getOuterBox(const StatePtr center,
                                    double size) const = 0;
          private:
            uint W_;
            std::vector<uint> lie_split_idx;

            inline void checkSetup(){
                BOOST_ASSERT_MSG(lie_split_idx.size()==W_, "setupManifold() was not called!");
            }

          protected:
            struct PrivilegedCoordinate{
                virtual std::string getName() const{ return {"<UntitledCoordinate>"}; }
                virtual unsigned int getWeight() const = 0;
                virtual TangentVector getTangent(const StatePtr center) const = 0;
                /* see if you need the following ... */
                /* virtual double inBoxCoefficient(double t) = 0;
                   virtual double outBoxCoefficient(double t) = 0;
                   virtual double compute(const StateType& center,
                                          const StateType& target) = 0;*/
            };

            void setupManifold(){
                // Initialize splitting sequence;
                W_ = 0;
                uint ci=0;
                for(auto& c: coordinates){
                  W_+=c->getWeight();
                  for(uint i=0;i<c->getWeight();i++){
                      lie_split_idx.push_back(ci);
                  }
                  ci++;
                }
                std::cout << "Splitting sequence initialized to: [";
                for (auto& i:lie_split_idx)
                  std::cout << i << ", ";
                std::cout << "\b\b];" << std::endl;
            }

        };

        class ReedsSheppManifold : public SubRiemannianManifold<ReedsSheppStateSpace>
        {

        private:
            typedef SubRiemannianManifold<ReedsSheppStateSpace> Base ;
            struct LongitudinalCoordinate : public Base::PrivilegedCoordinate{
                unsigned int getWeight() const
                {
                    return 1;
                }

                std::string getName() const
                {
                    return {"longitudinal"};
                }

                TangentVector getTangent(const StatePtr center) const
                {
                    auto center_ = center->as<StateType>();
                    return {cos(center_->getYaw()), sin(center_->getYaw()), 0.0};
                }
            };

            struct HeadingCoordinate : public Base::PrivilegedCoordinate{
                unsigned int getWeight() const
                {
                    return 1;
                }
                virtual std::string getName() const
                {
                    return {"heading"};
                }

                TangentVector getTangent(const StatePtr center) const
                {
                    return {0.0, 0.0, 1.0};
                }
            };

            struct LateralCoordinate : public Base::PrivilegedCoordinate{
                unsigned int getWeight() const
                {
                    return 2;
                }

                std::string getName() const
                {
                    return {"lateral"};
                }

                TangentVector getTangent(const StatePtr center) const
                {
                    auto center_ = center->as<StateType>();
                    return {-sin(center_->getYaw()), cos(center_->getYaw()), 0.0};
                }
            };

            static std::vector<double> stateToVec(const StatePtr x){
                return {x->as<StateType>()->getX(),
                        x->as<StateType>()->getY(),
                        x->as<StateType>()->getYaw()};
            }

            /*std::vector<double> getInBoxHalfSides(double T) const{
              double latsize = 4*rho_*(1-cos(T/rho_));
              double frontsize = T*(sqrt(3/2)-1);
              return {frontsize, latsize, T/rho_};
            }*/

            LateralCoordinate l_;
            LongitudinalCoordinate f_;
            HeadingCoordinate h_;

        public:

            ReedsSheppManifold(double rho_):
                SubRiemannianManifold(rho_)//, R_(rho_)
            {
                coordinates.push_back(&f_);
                coordinates.push_back(&h_);
                coordinates.push_back(&l_);

                this->setupManifold();
            }


            bool inPositiveHalfspace(const StatePtr state,
                                     const StatePtr pivot,
                                     const TangentVector normal) const{
                double dx = state->as<StateType>()->getX()-
                            pivot->as<StateType>()->getX();
                double dy = state->as<StateType>()->getY()-
                            pivot->as<StateType>()->getY();
                double dh = state->as<StateType>()->getYaw()-
                            pivot->as<StateType>()->getYaw();
                return dx*normal[0]+dy*normal[1]+dh*normal[2]>0;
            }

            class ReedsSheppOuterBox : public Base::OuterBox
            {
              private:
                const ReedsSheppManifold& M;
                double rho_;

                std::vector<double> getOutBoxHalfSides(double T) const{
                  double latsize;
                  if(T<rho_*M_PI/2){
                    latsize = rho_*(1-cos(T/rho_));
                  }else{
                    latsize = T+rho_*(1-M_PI/2); // R+T-R*PI/2
                  }
                  return {T, latsize, fmin(T/rho_, M_PI)};
                }

              public:
                ReedsSheppOuterBox(const StatePtr center_, double size_,
                                   const ReedsSheppManifold& manifold):
                    Base::OuterBox(center_, size_), M(manifold), rho_(M.rho_){}

                ~ReedsSheppOuterBox() = default;

                bool intersectsHyperplane(const StatePtr state,
                                      const std::vector<double>& normal) const
                {
                    auto center = ReedsSheppManifold::stateToVec(this->center);
                    auto pivot = ReedsSheppManifold::stateToVec(state);
                    auto halfsizes = getOutBoxHalfSides(size);

                    if(normal[0]==0 && normal[1]==0)
                        return fabs(center[2]-pivot[2])<halfsizes[2];

                    BOOST_ASSERT_MSG(normal[2]==0, "Unexpected normal vector for splitting plane");
                    // get the four points of the rectangle
                    TangentVector F(M.f_.getTangent(this->center));
                    ompl::utils::vec_scalar_multiply_inplace(F, halfsizes[0]);
                    TangentVector L(M.l_.getTangent(this->center));
                    ompl::utils::vec_scalar_multiply_inplace(L, halfsizes[1]);

                    double last_diff_n = 0;
                    for(int i=0;i<4;i++)
                    {
                        int a =  2*((i >> 0) & 1)-1;
                        int b = 2*((i >> 1) & 1)-1;
                        auto DF = ompl::utils::vec_scalar_multiply(F, a);
                        auto DL = ompl::utils::vec_scalar_multiply(L, b);
                        auto vertex = ompl::utils::vec_sum(ompl::utils::vec_sum(center, DF), DL);
                        double diff_n = ompl::utils::dot_product(ompl::utils::vec_diff(vertex, pivot), normal, 2);
                        if(last_diff_n*diff_n<0){
                            return true;
                        }

                        last_diff_n = diff_n;
                    }
                    return false;
                }
            };

            Base::OuterBoxPtr getOuterBox(const StatePtr center, double size) const {
                return OuterBoxPtr(new ReedsSheppOuterBox(center, size, *this));
            }
        };
    } // end of base namespace
} // end of ompl namespace

#endif
