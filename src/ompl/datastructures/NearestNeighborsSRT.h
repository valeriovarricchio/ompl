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

#include <set>
#include <list>

namespace ompl
{

// _MP meant to be a pointer to Motion (subclass of RRTStar)
// ManifoldType meant to be derived from SubRiemannianManifold
    template <typename _MP, typename ManifoldType>
    class NearestNeighborsSRT: public ompl::NearestNeighbors<_MP>
    {
        const ManifoldType& M;

        typedef typename ompl::NearestNeighbors<_MP> Base;

        class Node;
        typedef std::shared_ptr<Node> NodePtr;

        // Helper for KDTree structure
        struct Node{
            _MP motion;
            std::vector<double> normal;
            uint depth;
            std::array<NodePtr, 2> children;
            Node* parent;
            uint side;

            Node(const _MP& motion_,
                 const std::vector<double>& normal_):
                 motion(motion_), normal(normal_), depth(0)
            { }

            inline bool isLeaf(){
                return !(hasChild(0) || hasChild(1));
            }

            inline bool hasChild(uint which){
                return children[which].get();
            }

            void addChild(NodePtr c, bool side){
                children[side] = c;
                c->depth = depth+1;
                c->parent = this;
                c->side = side;
            }
        };

        // Helpers for a Bounded Priority Queue (BPQ)
        struct Qelem {
            _MP motion;
            double score;
            bool operator<(const Qelem& other) const{
                return score < other.score;
            }
        };

        struct BPQ {
            size_t kLim;
            double rLim;

            std::set<Qelem> Q;

            BPQ():
                kLim(std::numeric_limits<size_t>::infinity()),
                rLim(std::numeric_limits<double>::infinity())
            {}

            void setRlim(double rl_){
                rLim = rl_;
            }

            void setKlim(size_t kl_){
                kLim = kl_;
            }

            void insert(const Qelem& q){
                Q.insert(q);
                while(Q.rbegin()->score > rLim || Q.size() > kLim){
                    if(!Q.empty())
                        Q.erase(--Q.end());
                }
            }

            double getUpperBound(){
                if(rLim<std::numeric_limits<double>::infinity())
                    return rLim;

                if(Q.size()<kLim)
                    return std::numeric_limits<double>::infinity();

                return Q.rbegin()->score;
            }
        };

        NodePtr root;
        size_t size_;

    public:
        NearestNeighborsSRT(const ManifoldType& M_):
          ompl::NearestNeighbors<_MP>(),
          M(M_), root(NULL), size_(0)
        {
        }

        ~NearestNeighborsSRT() = default; // TODO is this ok?

        virtual void setDistanceFunction(const typename Base::DistanceFunction&)
        {
            OMPL_ERROR("The SRT algorithm assumes the distance function "
                       "associated to the subriemannian geometry.");
        }

        bool reportsSortedResults() const
        {
            return 1;
        }

        void clear(){
            root = NULL; // NOTE: this should destroy the entire tree, TODO check
        }

        void add(const _MP &data){
            if(!root.get()){
                root = NodePtr(
                   new Node(data, M.getSplittingNormal(data->state, 0))
                              );
                size_=1;
                return;
            }

            add(data, root);
        }

        void add(const _MP &data, const NodePtr top){   // TODO should be private
            int side = M.inPositiveHalfspace(data->state,
                                  top->motion->state,
                                  top->normal);

            if(top->hasChild(side)){
                add(data, top->children[side]);
            }else{
                NodePtr newnode(new Node(data,
                     M.getSplittingNormal(data->state, top->depth+1)));
                top->addChild(newnode, side);
                size_++;
            }
        }

        bool remove(const _MP &data){
            if(!root.get()){
                return false;
            }

            return remove(data, root);
        }

        bool remove(const _MP &data, NodePtr top){
            if(top->motion == data){ // TODO check if pointer comparison is right
                                     // or shall we compare states instead?
                if(!top->isLeaf()){
                    OMPL_ERROR("Cannot remove a non-leaf node!"); // TODO implement removal of non-leaves
                    return 0;
                }

                top->parent->children[top->side] = NULL; // this should effectively delete the node
                                                         // since no other references should exist
                                                         // TODO test
                if(size_ >0) // TODO maybe unnecessary (else unreachable)
                   size_--;
                return 1;
            }

            int side = M.inPositiveHalfspace(data->state,
                                  top->motion->state,
                                  top->normal);

            if(top->hasChild(side))
                return remove(data, top->children[side]);

            OMPL_ERROR("Asked to remove a node that was not found!");
            return 0;
        }

        _MP nearest(const _MP &data) const{
            BPQ Q;
            Q.setKlim(1);
            query(data, root, Q);
            return Q.Q.begin()->motion;
        }

        void nearestK(const _MP &data, std::size_t k, std::vector<_MP> &out) const{
            BPQ Q;
            Q.setKlim(k);
            query(data, root, Q);
            transcribeToOut(Q,out);
        }

        void nearestR(const _MP &data, double radius, std::vector<_MP> &out) const{
            BPQ Q;
            Q.setRlim(radius);
            query(data, root, Q);
            transcribeToOut(Q,out);
        }

        std::size_t size() const{
            return size_;
        }

        void list(std::vector<_MP> &data) const {
            data.clear();
            data.reserve(size_);
            listRecursion(root, data);
        }

    private:
        void listRecursion(NodePtr top, std::vector<_MP>& out) const{
            out.push_back(top->motion);
            if(top->hasChild(0))
                listRecursion(top->children[0], out);
            if(top->hasChild(1))
                listRecursion(top->children[1], out);
        }

        void query(const _MP &data, const NodePtr top, BPQ& Q) const{

            bool side = M.inPositiveHalfspace(data->state,
                                  top->motion->state,
                                  top->normal);

            if(top->hasChild(side)){
                query(data, top->children[side], Q);
            }

            Q.insert({top->motion, M.distance(data->state, top->motion->state)});

            if(top->hasChild(1-side)){
                auto oBox = M.getOuterBox(data->state, Q.getUpperBound());
                if(oBox->intersectsHyperplane(top->motion->state, top->normal))
                    query(data, top->children[1-side], Q);
            }
        }

        void transcribeToOut(BPQ Q, std::vector<_MP>& out) const{
            out.reserve(Q.Q.size());
            for(auto& a:Q.Q)
                out.push_back(a.motion);
        }

    };
}

#endif
