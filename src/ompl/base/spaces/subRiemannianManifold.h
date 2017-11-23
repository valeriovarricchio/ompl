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

/* Author: Valerio Varricchio */

#include <type_traits>

#ifndef OMPL_BASE_SPACES_SUB_RIEMANNIAN_MANIFOLD_
#define OMPL_BASE_SPACES_SUB_RIEMANNIAN_MANIFOLD_

namespace ompl {
namespace base {

/* _T must be a StateSpace type with a metric defined */
// TODO can we concept-check at compile time?
template <class _T>
class SubRiemannianManifold : public _T // TODO move out into SubRiemannianManifold.h
{

  using _T::_T; // inherit base ctors

  protected:
    struct PrivilegedCoordinate;

  public:

    typedef std::vector<double> TangentVector;
    std::vector<PrivilegedCoordinate*> coordinates;

    virtual ~SubRiemannianManifold(){} // TODO check this out

    virtual bool inPositiveHalfspace(const ompl::base::State* state,
                                     const ompl::base::State* pivot,
                                     const TangentVector normal) const=0;

    virtual TangentVector getSplittingNormal(const ompl::base::State* state,
                                             uint depth){
        checkSetup();
        return coordinates[lie_split_idx[depth%W_]]->getTangent(state);
    }

    class OuterBox{
      protected:
        double size;
        const ompl::base::State* center;
      public:
        OuterBox(const ompl::base::State* center_, double size_):
            size(size_), center(center_) // TODO note only copying pointer, relying on its validity!
        {}

        virtual bool intersectsHyperplane(const ompl::base::State* state,
                              const std::vector<double>& normal) const=0;

        virtual ~OuterBox() = default;
    };

    typedef std::shared_ptr<OuterBox> OuterBoxPtr;

    virtual OuterBoxPtr getOuterBox(const ompl::base::State* center, double size) const=0;
  private:
    uint W_;
    std::vector<uint> lie_split_idx;

    inline void checkSetup(){
        BOOST_ASSERT_MSG(lie_split_idx.size()==W_,
                         "setupManifold() was not called!");
    }

  protected:
    struct PrivilegedCoordinate{
        virtual std::string getName() const{
            return {"<UntitledCoordinate>"};
        }
        virtual unsigned int getWeight() const = 0;
        virtual TangentVector getTangent(const ompl::base::State* center) const=0;
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

} // end namespace base
} // end namespace ompl

#endif

