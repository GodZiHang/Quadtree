#ifndef KDTREE_H
#define KDTREE_H

#include <cassert>
#include <algorithm>
#include <array>
#include <memory>
#include <type_traits>
#include <vector>
#include "box.h"

/**
 * This file is a template implementation of a two-dimensional kd-tree,
 * which differs from a regular kd-tree in the following ways:
 *   - The node it saves is not a point, but a rectangular area
 *   - Its regional division direction is not alternating, and during
 *   initialization, a priority axis will be passed in:
 *     + When depth<=maxDepth/2, divide according to the priority axis direction
 *     + When depth>maxDepth/2, divide in another direction
 *   - region divides by the midpoint of the division direction, rather
 *    than taking the coordinates of the intermediate elements for division
 *-----------------------------------------------------------------------------
 * 该文件是一个二维kd-tree的模板实现，它与普通的kd-tree有以下不同：
 *   - 它保存的节点不是一个点，而是一个矩形区域
 *   - 它的区域划分方向并非交替进行，在初始化时，会传入一个优先轴：
 *     + 在 depth <= maxDepth/2 时，按优先轴方向进行划分
 *     + 在 depth >  maxDepth/2 时，按另一个方向进行划分
 *   - 每次进行区域划分时，它按划分方向的中点进行二分，而不是取中间
 *   元素的坐标进行划分
 */
namespace kdTree
{
using namespace quadtree;

template<typename T, typename GetBox, typename Equal = std::equal_to<T>, typename Float = float>
class kdTree_2
{
    static_assert(std::is_convertible_v<std::invoke_result_t<GetBox, const T&>, Box<Float>>,
        "GetBox must be a callable of signature Box<Float>(const T&)");
    static_assert(std::is_convertible_v<std::invoke_result_t<Equal, const T&, const T&>, bool>,
        "Equal must be a callable of signature bool(const T&, const T&)");
    static_assert(std::is_arithmetic_v<Float>);

public:
    enum priorAxis {xFirst, yFirst};
    kdTree_2(const Box<Float>& box, const GetBox& getBox = GetBox(),
             priorAxis prior = xFirst, const Equal& equal = Equal()) :
             mBox(box), mRoot(std::make_unique<Node>()), mGetBox(getBox), mEqual(equal), mPrior(prior)
         {

         }

    void add(const T& value)
    {
        add(mRoot.get(), 0, mBox, value);
    }

    void remove(const T& value)
    {
        remove(mRoot.get(), mBox, value);
    }

    std::vector<T> query(const Box<Float>& box) const
    {
        auto values = std::vector<T>();
        query(mRoot.get(), mBox, box, values);
        return values;
    }

    void effectiveQuery(const Box<Float>& box, T* valuesArr, int& arrSize, int& maxSize)
    {
        effectiveQuery(mRoot.get(), mBox, box, valuesArr, arrSize, maxSize);
    }

    Box<Float> getBox() const
    {
        return mBox;
    }

private:
    static constexpr auto Threshold = std::size_t(16);
    static constexpr auto MaxDepth = std::size_t(20);
    enum Orient { Smaller = 0, Bigger = 1, Invalid = -1 };
    struct Node
    {
        std::array<std::unique_ptr<Node>, 2> children;
        std::vector<T> values;
        int axis = -1; // parting axis, -1 for leaf, 0 for x, 1 for y
    };

    Box<Float> mBox;
    std::unique_ptr<Node> mRoot;
    GetBox mGetBox;
    Equal mEqual;
    priorAxis mPrior;

    bool isLeaf(const Node* node) const
    {
        return !static_cast<bool>(node->children[0]);
    }

    Box<Float> computeBox(const Box<Float>& box, Orient i, int axis) const
    {
        if(axis == 0 && i == Orient::Smaller) return Box<Float>(box.left, box.top, box.width / 2, box.height);
        if(axis == 0 && i == Orient::Bigger) return Box<Float>(box.left + box.width / 2, box.top, box.width / 2, box.height);
        if(axis == 1 && i == Orient::Smaller) return Box<Float>(box.left, box.top, box.width, box.height / 2);
        if(axis == 1 && i == Orient::Bigger) return Box<Float>(box.left, box.top + box.height / 2, box.width, box.height / 2);

        // Should never happen
        assert(false && "Invaild axis or i");
        return Box<Float>();
    }

    Orient getOrient(const Box<Float>& nodeBox, int axis, const Box<Float>& valueBox) const
    {
        assert(axis != -1);
        auto center = nodeBox.getCenter();
        if(axis == 0)
        {
            if(valueBox.getRight() < center.x) return Orient::Smaller;
            else if(valueBox.left >= center.x) return Orient::Bigger;
            else return Orient::Invalid;
        }
        else if(axis == 1)
        {
            if(valueBox.getBottom() < center.y) return Orient::Smaller;
            else if(valueBox.top >= center.y) return Orient::Bigger;
            else return Orient::Invalid;
        }
        else
        {
            assert(false && "Invalid axis");
            return Orient::Invalid;
        }
    }
    
    void add(Node* node, std::size_t depth, const Box<Float>& box, const T& value)
    {
        assert(node != nullptr);
        assert(box.contains(mGetBox(value)));
        if(isLeaf(node))
        {
            node->values.push_back(value);
            if(node->values.size() > Threshold && depth < MaxDepth)
            {
                split(node, depth, box);
            }
        }
        else
        {
            auto i = getOrient(box, node->axis, mGetBox(value));
            if(i == Orient::Smaller)
            {
                add(node->children[0].get(), depth + 1, computeBox(box, Orient::Smaller, node->axis), value);
            }
            else if(i == Orient::Bigger)
            {
                add(node->children[1].get(), depth + 1, computeBox(box, Orient::Bigger, node->axis), value);
            }
            else
            {
                node->values.push_back(value);
            }
        }
    }

    void split(Node* node, std::size_t depth, const Box<Float>& box)
    {
        assert(node != nullptr);
        assert(isLeaf(node) && "Only leaf nodes can be split");
        // create children
        node->children[0] = std::make_unique<Node>();
        node->children[1] = std::make_unique<Node>();
        if(depth <= MaxDepth / 2) node->axis = mPrior;
        else node->axis = 1 ^ mPrior; // Used to reduce branches, only applicable to the current case
        // assign values to children
        auto newValues = std::vector<T>();
        for(const auto& value : node->values) 
        {
            auto i = getOrient(box, node->axis, mGetBox(value));
            if(i != Orient::Invalid)
                node->children[static_cast<std::size_t>(i)]->values.push_back(value);
            else
                newValues.push_back(value);
        }
        node->values = std::move(newValues);
    }

    bool remove(Node* node, const Box<Float>& box, const T& value)
    {
        assert(node != nullptr);
        assert(box.contains(mGetBox(value)));
        if(isLeaf(node))
        {
            removeValue(node, value);
            return true;
        }
        else 
        {
            auto i = getOrient(box, node->axis, mGetBox(value));
            if(i == 0) {
                if(remove(node->children[static_cast<std::size_t>(i)].get(), computeBox(box, Orient::Smaller, node->axis), value))
                    return tryMerge(node);
            }
            else if(i == 1) {
                if(remove(node->children[static_cast<std::size_t>(i)].get(), computeBox(box, Orient::Bigger, node->axis), value))
                    return tryMerge(node);
            }
            else removeValue(node, value);
        }
        return false;
    }

    void removeValue(Node* node, const T& value)
    {
        assert(node != nullptr);
        auto it = std::find_if(node->values.begin(), node->values.end(), 
            [&](const auto& v) { return mEqual(v, value); });
        if(it != node->values.end())
        {
            node->values.erase(it);
        }
    }

    bool tryMerge(Node* node)
    {
        assert(node != nullptr);
        assert(!isLeaf(node));
        if(isLeaf(node->children[0].get()) && isLeaf(node->children[1].get()))
        {
            std::size_t nValues = node->values.size() + node->children[0]->values.size() + node->children[1]->values.size();
            if(nValues <= Threshold)
            {
                node->values.reserve(nValues);
                node->values.insert(node->values.end(), node->children[0]->values.begin(), node->children[0]->values.end());
                node->values.insert(node->values.end(), node->children[1]->values.begin(), node->children[1]->values.end());
                node->children[0].reset();
                node->children[1].reset();
                return true;
            }
            return false;
        }
        return false;
    }

    void query(Node* node, const Box<Float>& box, const Box<Float> queryBox, std::vector<T>& values) const
    {
        assert(node != nullptr);
        assert(queryBox.intersects(box));
        if(queryBox.contains(box)) 
        {
            values.insert(values.end(), node->values.begin(), node->values.end());
        }
        else
        {
            for(const auto& value : node->values)
            {
                if(queryBox.intersects(mGetBox(value)))
                    values.push_back(value);
            }
        }
        
        if(!isLeaf(node))
        {
            auto childBox = computeBox(box, Orient::Smaller, node->axis);
            if(queryBox.intersects(childBox))
                query(node->children[0].get(), childBox, queryBox, values);
            childBox = computeBox(box, Orient::Bigger, node->axis);
            if(queryBox.intersects(childBox))
                query(node->children[1].get(), childBox, queryBox, values);
        }
    }

    void effectiveQuery(Node* node, const Box<Float>& box, const Box<Float>& queryBox, T* valuesArr, int& arrSize, int& maxSize)
    {
        if(queryBox.contains(box)) {
            effectiveQuery_addAllElements(node, valuesArr, arrSize, maxSize);
            return;
        } else {
            for(const auto& value: node->values)
            {
                if(queryBox.intersects(mGetBox(value)))
                    valuesArr[arrSize++] = value;
            }
        }

        if(!isLeaf(node))
        {
            auto childBox = computeBox(box, Smaller, node->axis);
            if (queryBox.intersects(childBox)) {
                effectiveQuery(node->children[0].get(), childBox, queryBox, valuesArr, arrSize, maxSize);
            }
            childBox = computeBox(box, Bigger, node->axis);
            if (queryBox.intersects(childBox)) {
                effectiveQuery(node->children[1].get(), childBox, queryBox, valuesArr, arrSize, maxSize);
            }
        }
    }

    void effectiveQuery_addAllElements(Node* node, T* valuesArr, int& arrSize, int& maxSize)
    {
        for (const auto& value : node->values)
            valuesArr[arrSize++] = value;
        if(!isLeaf(node))
        {
            for (auto i = std::size_t(0); i < node->children.size(); ++i)
            {
                effectiveQuery_addAllElements(node->children[i].get(), valuesArr, arrSize, maxSize);
            }
        }
    }
};
}

#endif // KDTREE_H

