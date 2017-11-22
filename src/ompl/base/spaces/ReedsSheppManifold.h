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
#include "ompl/base/spaces/subRiemannianManifold.h"
#include <boost/math/constants/constants.hpp>

namespace ompl {
namespace base {

class ReedsSheppManifold :
        public SubRiemannianManifold<ReedsSheppStateSpace> {

private:
    typedef SubRiemannianManifold<ReedsSheppStateSpace> Base ;
    struct LongitudinalCoordinate : public Base::PrivilegedCoordinate{
        unsigned int getWeight() const {
            return 1;
        }

        std::string getName() const {
            return {"longitudinal"};
        }

        TangentVector getTangent(const StatePtr center) const {
            auto x = center->as<StateType>();
            return {cos(x->getYaw()), sin(x->getYaw()), 0.0};
        }
    };

    struct HeadingCoordinate : public Base::PrivilegedCoordinate{
        unsigned int getWeight() const {
            return 1;
        }
        virtual std::string getName() const {
            return {"heading"};
        }

        TangentVector getTangent(const StatePtr center) const {
            return {0.0, 0.0, 1.0};
        }
    };

    struct LateralCoordinate : public Base::PrivilegedCoordinate{
        unsigned int getWeight() const {
            return 2;
        }

        std::string getName() const {
            return {"lateral"};
        }

        TangentVector getTangent(const StatePtr center) const {
            auto x = center->as<StateType>();
            return {-sin(x->getYaw()), cos(x->getYaw()), 0.0};
        }
    };

    LateralCoordinate l_;
    LongitudinalCoordinate f_;
    HeadingCoordinate h_;

public:

    ReedsSheppManifold(double rho_): SubRiemannianManifold(rho_)
    {
        coordinates.push_back(&f_);
        coordinates.push_back(&h_);
        coordinates.push_back(&l_);

        this->setupManifold();
    }

    bool inPositiveHalfspace(const StatePtr state, const StatePtr pivot,
                             const TangentVector normal) const
    {
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
        typedef std::array<double, 2> vec2d;

        const ReedsSheppManifold& M;
        double rho_;

        std::vector<double> getOutBoxHalfSides(double T) const
        {
          double latsize;
          if(T<rho_*M_PI/2){
            latsize = rho_*(1-cos(T/rho_));
          }else{
            latsize = T+rho_*(1-M_PI/2);
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
            if(std::isinf(size))
                return true;

            std::vector<double> halfsizes(getOutBoxHalfSides(size));

            if(normal[0]==0 && normal[1]==0)
                return fabs(center->as<StateType>()->getYaw()
                            -state->as<StateType>()->getYaw())
                        < halfsizes[2];

            BOOST_ASSERT_MSG(normal[2]==0,
                    "Unexpected normal vector for splitting plane");

            vec2d rcenter({center->as<StateType>()->getX()
                           -state->as<StateType>()->getX(),
                           center->as<StateType>()->getY()
                           -state->as<StateType>()->getY()});

            TangentVector f(M.f_.getTangent(this->center));
            TangentVector l(M.l_.getTangent(this->center));
            vec2d F({f[0]*halfsizes[0], f[1]*halfsizes[0]});
            vec2d L({l[0]*halfsizes[1], l[1]*halfsizes[1]});

            double last_dotp = 0;
            int a,b;
            for(int i=0;i<4;i++)
            {
                a = ((i & 1) << 1)-1;
                b = (i & 2)-1;
                vec2d rvertex({rcenter[0]+F[0]*a+L[0]*b,
                             rcenter[1]+F[1]*a+L[1]*b});
                double dotp = rvertex[0]*normal[0]+rvertex[1]*normal[1];

                if(last_dotp*dotp<0)
                    return true;

                last_dotp = dotp;
            }
            return false;
        }
    };

    OuterBoxPtr getOuterBox(const StatePtr x0, double size) const {
        return OuterBoxPtr(new ReedsSheppOuterBox(x0, size, *this));
    }
};

} // end of base namespace
} // end of ompl namespace

#endif
