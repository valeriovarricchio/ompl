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

#include <easy/profiler.h>

namespace ompl {
namespace base {

class ReedsSheppManifold :
        public SubRiemannianManifold<ReedsSheppStateSpace> {

private:
    typedef std::array<double, 2> vec2d;
    typedef SubRiemannianManifold<ReedsSheppStateSpace> Base ;
    struct LongitudinalCoordinate : public Base::PrivilegedCoordinate{
        unsigned int getWeight() const {
            return 1;
        }

        std::string getName() const {
            return {"longitudinal"};
        }

        TangentVector getTangent(const ompl::base::State* center) const {
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

        TangentVector getTangent(const ompl::base::State* center) const {
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

        TangentVector getTangent(const ompl::base::State* center) const {
            auto x = center->as<StateType>();
            return {-sin(x->getYaw()), cos(x->getYaw()), 0.0};
        }
    };

    template<class at, class bt>
    static double dotXY(const at& A, const bt& B) {
        return A[0]*B[0]+A[1]*B[1];
    }

    static inline double mod2pi(const double& v) {
        constexpr double pi = boost::math::constants::pi<double>();
        constexpr double twopi = 2.0*pi;
        double rsp = fmod(v, twopi);
        if (rsp < -pi)
            rsp += twopi;
        else if (rsp > pi)
            rsp -= twopi;
        return rsp;
    }

    static inline vec2d diffXY(const ompl::base::State* A,
                               const ompl::base::State* B) {
        return {A->as<StateType>()->getX()-B->as<StateType>()->getX(),
                A->as<StateType>()->getY()-B->as<StateType>()->getY()};
    }

    double inline euclideanXY(const ompl::base::State* A,
                      const ompl::base::State* B) const {
        vec2d dxy(diffXY(A,B));
        return sqrt(dxy[0]*dxy[0]+dxy[1]*dxy[1]);
    }

    double df(const ompl::base::State* A, const ompl::base::State* B) const{
        vec2d dxy(diffXY(A,B));
        return std::max(fabs(dotXY(f_.getTangent(A), dxy)),
                        fabs(dotXY(f_.getTangent(B), dxy)));
    }

    double dl(const ompl::base::State* A, const ompl::base::State* B) const{
        vec2d dxy(diffXY(A,B));
        return std::max(fabs(dotXY(l_.getTangent(A), dxy)),
                        fabs(dotXY(l_.getTangent(B), dxy)));
    }

    double dh(const ompl::base::State* A, const ompl::base::State* B) const{
        return std::min(fabs(mod2pi(A->as<StateType>()->getYaw()-
                                    B->as<StateType>()->getYaw())), M_PI);
    }

    void reduceRange(const srt::Node& n, const double yaw, vec2d& thetaMinMax) const{
        if(n.normal[2]==1){
            double nodeyaw = n.state->as<StateType>()->getYaw();
            bool side = yaw>=nodeyaw;
            bool comp = (thetaMinMax[1-side]-nodeyaw)>0;
            if(comp!=side){
                thetaMinMax[1-side] = nodeyaw;
            }
        }
    }

    double thetaRange(const srt::Bucket& bucket, const double yaw) const{
        vec2d thetaMinMax({-M_PI, M_PI});
        if(useCache_){ // root to leaf
            auto cur = bucket.nodes.rbegin(); // root
            while(cur != bucket.nodes.rend() && cur->cached) cur++;
            if(cur!=bucket.nodes.rbegin()){
                thetaMinMax=*static_cast<vec2d*>(std::prev(cur)->cached);
                reduceRange(*std::prev(cur), yaw, thetaMinMax);
            }
            while(cur != bucket.nodes.rend()){
                cur->cached = new vec2d(thetaMinMax);
                reduceRange(*cur, yaw, thetaMinMax);
                cur++;
            }
        }else{ // leaf to root
            for (const auto& node:bucket.nodes)
                reduceRange(node, yaw, thetaMinMax);
        }
        assert(thetaMinMax[1]-thetaMinMax[0]>0 && "WRONG THETA RANGE!");
        return thetaMinMax[1]-thetaMinMax[0];
    }


    LateralCoordinate l_;
    LongitudinalCoordinate f_;
    HeadingCoordinate h_;
    const bool useTransition_; // TODO can these be known at compile time? (constexpr)
    const bool useLowerBound_;
    const bool useCache_;
    const double transThreshold_;
public:

    ReedsSheppManifold(double rho_,
                       bool useTransition=true,
                       bool useLowerBound=true,
                       bool useCache=true,
                       double transThreshold = M_PI):
        SubRiemannianManifold(rho_),
        useTransition_(useTransition),
        useLowerBound_(useLowerBound),
        useCache_(useCache),
        transThreshold_(transThreshold)
    {
        coordinates.push_back(&f_);
        coordinates.push_back(&h_);
        coordinates.push_back(&l_);

        this->setupManifold();
    }

    bool inPositiveHalfspace(const ompl::base::State* state,
                             const ompl::base::State* pivot,
                             const TangentVector& normal) const
    {
        EASY_BLOCK("srt: inPositiveHalfspace");
        auto s = state->as<StateType>();
        auto p = pivot->as<StateType>();
        double dx = s->getX()-p->getX();
        double dy = s->getY()-p->getY();
        double dh = s->getYaw()-p->getYaw();
        return dx*normal[0]+dy*normal[1]+dh*normal[2]>0;
    }

    inline bool hasLowerBound() const{
        return useLowerBound_;
    }

    bool needsGhostPoints(const ompl::base::State* q, double dmax) const override {
        int m = (q->as<StateType>()->getYaw() >0)? -1 : 1;
        double dthreshold = (M_PI + m*q->as<StateType>()->getYaw())*rho_;
        return dmax>dthreshold;
    }

    double lowerBound(const ompl::base::State* A,
                      const ompl::base::State* B) const {
        EASY_BLOCK("srt: lowerBound");
        double tf = df(A,B);
        double th = dh(A,B)*rho_;
        double DL = dl(A,B);
        double tl = (DL>rho_)?
                      (DL-rho_)+M_PI*rho_/2:
                      acos(1-DL/rho_)*rho_;
        double tbox = std::max(std::max(tf, th), tl);
        return std::max(tbox, euclideanXY(A,B));
        //return tbox; // only use box bound
    }

    /* This implements the transition */
    TangentVector getSplittingNormal(const ompl::base::State* state,
                                     const srt::Bucket& bucket) const {
        EASY_BLOCK("srt: getNormal");
        if(useTransition_){
            uint depth = bucket.nodes.size();
            double thrange = thetaRange(bucket, state->as<StateType>()->getYaw());

            if(thrange>transThreshold_){
                // "holonomic regime"
                //return coordinates[depth%3]->getTangent(state);
                TangentVector rsp({0.0,0.0,0.0});
                rsp[depth%3] = 1.0;
                return rsp;
            }
        }
        // "non-holonomic regime"
        return Base::getSplittingNormal(state, bucket);
    }

    void freeCacheData(void* cacheData) const{
        delete static_cast<vec2d*>(cacheData);
    }

    class ReedsSheppOuterBox : public Base::OuterBox
    {
      private:

        const ReedsSheppManifold& M;
        double rho_;

        std::vector<double> getOutBoxHalfSides(double T) const
        {
          double latsize;
          if(T<rho_*M_PI/2){             // TODO use boost pi?
            latsize = rho_*(1-cos(T/rho_));
          }else{
            latsize = T+rho_*(1-M_PI/2);
          }
          return {T, latsize, fmin(T/rho_, M_PI)};
        }

      public:
        ReedsSheppOuterBox(const ompl::base::State* center_, double size_,
                           const ReedsSheppManifold& manifold):
            Base::OuterBox(center_, size_), M(manifold), rho_(M.rho_){}

        ~ReedsSheppOuterBox() = default;

        bool intersectsHyperplane(const ompl::base::State* state,
                                  const std::vector<double>& normal) const
        {
            EASY_BLOCK("srt: ballHyperplane");
            if(std::isinf(size))
                return true;

            if(size < fabs(dotXY(diffXY(center, state), normal)))
                return false;

            std::vector<double> halfsizes(getOutBoxHalfSides(size));

            if(normal[0]==0 && normal[1]==0)
                return fabs(center->as<StateType>()->getYaw()
                            -state->as<StateType>()->getYaw())
                        < halfsizes[2];

            BOOST_ASSERT_MSG(normal[2]==0,
                    "Unexpected normal vector for splitting plane");

            vec2d rcenter(diffXY(center, state));

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

    OuterBoxPtr getOuterBox(const ompl::base::State* x0, double size) const {
        return OuterBoxPtr(new ReedsSheppOuterBox(x0, size, *this));
    }
};

} // end of base namespace
} // end of ompl namespace

#endif
