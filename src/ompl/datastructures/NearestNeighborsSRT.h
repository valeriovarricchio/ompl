/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2008, Willow Garage, Inc.
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
*   * Neither the name of the Willow Garage nor the names of its
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

#ifndef OMPL_DATASTRUCTURES_NEAREST_NEIGHBORS_SRT_
#define OMPL_DATASTRUCTURES_NEAREST_NEIGHBORS_SRT_

#include <vector>
#include <functional>
#include <ompl/base/StateSpace.h>
#include <ompl/datastructures/NearestNeighbors.h>

namespace ompl
{

// _MP meant to be a pointer to Motion (subclass of RRTStar)
// ManifoldType meant to be derived from SubRiemannianManifold
    template <typename _MP, typename ManifoldType>
    class NearestNeighborsSRT: public ompl::NearestNeighbors<_MP>
    {
        const ManifoldType& M;

        typedef typename ompl::NearestNeighbors<_MP> Base;

        // Define private sub classes

    public:
        NearestNeighborsSRT(const ManifoldType& M_):
          ompl::NearestNeighbors<_MP>(),
          M(M_)
        {
        }

        ~NearestNeighborsSRT() = default;

        virtual void setDistanceFunction(const typename Base::DistanceFunction&)
        {
            OMPL_ERROR("The SRT algorithm assumes the distance function "
                       "associated to the subriemannian geometry.");
        }

        bool reportsSortedResults() const
        {
            return 1;
        }

        void clear(){ // TODO
        }

        void add(const _MP &data){
            std::cout << "add() called" << std::endl;
            // TODO
        }

        bool remove(const _MP &data){
            std::cout << "remove() called" << std::endl;
            OMPL_ERROR("remove() hasn't been implemented in STR.");
            return 0;
        }

        _MP nearest(const _MP &data) const{
            std::cout << "nearest() called" << std::endl;
            return {};
        }

        void nearestK(const _MP &data, std::size_t k, std::vector<_MP> &nbh) const{
            std::cout << "nearestK() called" << std::endl;
        }

        void nearestR(const _MP &data, double radius, std::vector<_MP> &nbh) const{
            std::cout << "nearestR() called" << std::endl;
            OMPL_ERROR("Range search hasn't been implemented in STR.");
        }

        std::size_t size() const{
            std::cout << "size() called" << std::endl;
            return 0;
        }

        /** \brief Get all the elements in the datastructure */
        void list(std::vector<_MP> &data) const {
            std::cout << "list() called" << std::endl;
        }
    };
}

#endif
