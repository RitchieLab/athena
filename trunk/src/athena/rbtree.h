/*
Copyright Marylyn Ritchie 2011

This file is part of ATHENA.

ATHENA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ATHENA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ATHENA.  If not, see <http://www.gnu.org/licenses/>.
*/
//---------------------------------------------------------------------------
#ifndef RBTREE_H
#define RBTREE_H

//---------------------------------------------------------------------------
// Template implementation
// Red-Black Tree with Doubly Linked List
//
// Author:  Arjan van den Boogaard (arjan@threelittlewitches.nl)
//          Three Little Witches (www.threelittlewitches.nl)
//          7 feb 2004
//
// Modifications:
// 			Eric Torstenson
//			Templated the comparison function
//
// This is free software and you may use it as you wish, at your own risk
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
#include <stdlib.h>
#include <iostream>
#include <assert.h>

namespace Utility {

template <class K>
struct rbLT {
	int operator()(const K l, const K r) const {
		if (l<r) return -1;
		if (l>r) return 1;
		return 0;
	}
};
	
//---------------------------------------------------------------------------
template <class K, class D, typename C = rbLT<K> > class RBTree;
template <class K, class D, typename C = rbLT<K> > class RBTreeNode;

//---------------------------------------------------------------------------
enum RBTreeColor { scRED, scBLACK, scDUP };



//---------------------------------------------------------------------------
template <class K, class D, typename C > class RBTreeNode
{
    public:
        // Returns the key. The key is read-only but can be changed by
        // using the ReplaceKey method in the storage object
        K GetKey()                          { return Key;  }

        // Get or set the user data
        D       GetData()                   { return Data; }
        void    SetData(const D &newdata)   { Data = newdata; }

        // Get the previous or next node in the storage. This ordered
        // traversal includes the duplicates. For only the unique nodes use
        // the GetNextUnique and GetPrevUnique methods in the storage object
        RBTreeNode<K,D,C> *  GetNext()      { return Next; }
        RBTreeNode<K,D,C> *  GetPrev()      { return Prev; }

    private:
        friend class RBTree<K,D,C>;

        RBTreeNode<K,D,C> * Next;
        RBTreeNode<K,D,C> * Prev;
        RBTreeNode<K,D,C> * Left;
        RBTreeNode<K,D,C> * Right;
        RBTreeNode<K,D,C> * Parent;
        RBTreeColor       Color;

        K   Key;
        D   Data;
};




//---------------------------------------------------------------------------
template <class K, class D, typename C > class RBTree
{
    public:
        RBTree();

		// Use this ONLY if you want to limit the size of the tree. It will always prefer 
		// CompareFn() > 0 and will delete nodes if necessary to maintain the size
		RBTree(int maxSize);

		RBTree(const RBTree& other);

        // Deleting a storage object clears all remaining nodes
        ~RBTree();

		RBTree& operator=(const RBTree& other);

        // Add and Insert are the same (for backwards compatibility :)
        RBTreeNode<K,D,C> * Insert(const K &key, const D &data);
        RBTreeNode<K,D,C> * Add(const K &key, const D &data);

        // Deletes the specified node from the tree. The pointer is no
        // longer valid after the delete
        void Delete(RBTreeNode<K,D,C> * node);

        // Find the first instance of the specified value
        RBTreeNode<K,D,C> * Find(const K &key);
		RBTreeNode<K,D,C> * FindNearestMin(const K &key);
		RBTreeNode<K,D,C> * FindNearestMax(const K &key);

        // Replaces the key (and thus change sorting order)
        void  ReplaceKey(RBTreeNode<K,D,C> * node, const K &newkey);

        // Clear wipes the complete tree
        void Clear();
        
        // Sets maximum size of tree
        void SetMaxSize(int max);

        // First and Last are the head and tail of the doubly linked list
        RBTreeNode<K,D,C> * GetFirst() const;
        RBTreeNode<K,D,C> * GetLast() const;

        // Returns the total amount of nodes in the storage
        int GetCount();

        // Traversal through all unique nodes
        RBTreeNode<K,D,C> * GetNextUnique(RBTreeNode<K,D,C> * node);
        RBTreeNode<K,D,C> * GetPrevUnique(RBTreeNode<K,D,C> * node);

    protected:
		// The compare function MUST return -1, 0, 1. The tree is built upon this and searches
		// must behave the same. 
 		C 	ComparFn;

        RBTreeNode<K,D,C> *  First;
        RBTreeNode<K,D,C> *  Last;
        RBTreeNode<K,D,C> *  Root;
        RBTreeNode<K,D,C> *  NIL;
        int                  mCount;
		int 				 maxSize;

        void RotateRight(RBTreeNode<K,D,C> * x);
        void RotateLeft(RBTreeNode<K,D,C> * x);
        void FixInsert(RBTreeNode<K,D,C> * x);
        void FixDelete(RBTreeNode<K,D,C> * x);
        void DoInsert(RBTreeNode<K,D,C> * x);
        void DoDelete(RBTreeNode<K,D,C> * x);

        void DLLDelete(RBTreeNode<K,D,C> * x);
        void DLLInsertBefore(RBTreeNode<K,D,C> * x, RBTreeNode<K,D,C> * y);
        void DLLAddAfter(RBTreeNode<K,D,C> * x, RBTreeNode<K,D,C> * y);

        RBTreeNode<K,D,C> *  TreeNext(RBTreeNode<K,D,C> * x);
        RBTreeNode<K,D,C> *  TreePrev(RBTreeNode<K,D,C> * x);
        RBTreeNode<K,D,C> *  FindMin(RBTreeNode<K,D,C> * x);
        RBTreeNode<K,D,C> *  FindMax(RBTreeNode<K,D,C> * x);

        void InitNode(RBTreeNode<K,D,C> * x);
        void ReplaceTreeNode(RBTreeNode<K,D,C> * x, RBTreeNode<K,D,C> * y);
};



//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
RBTree<K,D,C>::RBTree() : mCount(0), maxSize(0)
{
    NIL = new RBTreeNode<K,D,C>;
    InitNode(NIL);
    NIL->Color = scBLACK;

    Root = NIL;
    First = NULL;
    Last = NULL;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
RBTree<K,D,C>::RBTree(int maxSize) : mCount(0), maxSize(maxSize)
{
    NIL = new RBTreeNode<K,D,C>;
    InitNode(NIL);
    NIL->Color = scBLACK;

    Root = NIL;
    First = NULL;
    Last = NULL;
}

template <class K, class D, typename C>
inline
RBTree<K,D,C>::RBTree(const RBTree<K,D,C>& other) : mCount(other.mCount), maxSize(other.maxSize) {
	NIL = new RBTreeNode<K,D,C>;
    InitNode(NIL);
    NIL->Color = scBLACK;

    Root = NIL;
    First = NULL;
    Last = NULL;

	RBTreeNode<K,D,C> *node = other.GetFirst();
	while (node) {
		Insert(node->GetKey(), node->GetData());
		node = node->GetNext();
	}	
}

template <class K, class D, typename C>
inline
RBTree<K,D,C>& RBTree<K,D,C>::operator=(const RBTree<K,D,C>& other) {
    Clear();
    delete NIL;

	NIL = new RBTreeNode<K,D,C>;
    InitNode(NIL);
    NIL->Color = scBLACK;

    Root = NIL;
    First = NULL;
    Last = NULL;

	RBTreeNode<K,D,C> *node = other.GetFirst();
	while (node) {
		Insert(node->GetKey(), node->GetData());
		node = node->GetNext();
	}	
	return *this;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
RBTree<K,D,C>::~RBTree()
{
    Clear();
    delete NIL;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::Add(const K &key, const D &data)
{
    return Insert(key, data);
}


//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
void RBTree<K,D,C>::SetMaxSize(int max)
{
    maxSize = max;
}

// #include <iostream>
// using namespace std;

//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::Insert(const K &key, const D &data)
{
    RBTreeNode<K,D,C> *  node = NULL;
	/** 
	 * @todo On Insert where maxSize > 0, I'm using FindMin() to check whether or not
	 * I want to add the node. It might be better to just add it and then delete or
	 * to simply record the minimum during insertion/deletion
	 */
	if (maxSize == 0 || maxSize > mCount || key > FindMin(Root)->Key){
		node = new RBTreeNode<K,D,C>;
		InitNode(node);
	
		node->Key = key;
		node->Data = data;
		DoInsert(node);
	
		mCount++;

		while (maxSize > 0 && mCount > maxSize)
			Delete(FindMin(Root));
	  }
    return node;
}


//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
void RBTree<K,D,C>::Delete(RBTreeNode<K,D,C> * node)
{
    assert(node);

    DoDelete(node);

    mCount--;

    delete node;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::Find(const K &key)
{
    RBTreeNode<K,D,C> *  y = Root;
    int c;

    while (y != NIL)
    {
        c = ComparFn(key, y->Key);
        if (c == 0) return y;
        if (c < 0)  y = y->Left;
        else        y = y->Right;
    }

    return NULL;
}

template <class K, class D, typename C >
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::FindNearestMin(const K &key)
{
    RBTreeNode<K,D,C> *  y = Root;
	RBTreeNode<K,D,C> *  last = NULL;
    int c;

	while (y != NIL)
    {
       	c = ComparFn(key, y->Key);
        if (c == 0) return y;
        if (c < 0)   {
			last = y;
			y = y->Left; 	
		}
        else {
	        y = y->Right;
		}
    }

    return last;
}

#include <iostream> // remove
using namespace std; // remove

template <class K, class D, typename C >
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::FindNearestMax(const K &key)
{
    RBTreeNode<K,D,C> *  y = Root;
	RBTreeNode<K,D,C> *  last = NULL;
	
last = Root;	
    int c;
    // switched the last to be on the right side instead of the left side -- smd 5/22/09
    while (y != NIL)
    {
       	c = ComparFn(key, y->Key);
        if (c == 0) return y;
        if (c < 0)   {
			y = y->Left; 
		}
        else {
          last = y;
	        y = y->Right;
		}
    }
    return last;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C >
inline
void RBTree<K,D,C>::ReplaceKey(RBTreeNode<K,D,C> * node, const K &newkey)
{
    assert(node);

    DoDelete(node);
    InitNode(node);

    node->Key = newkey;

    DoInsert(node);
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::Clear()
{
    RBTreeNode<K,D,C> *  x = First;
    RBTreeNode<K,D,C> *  y;

    while (x)
    {
        y = x;
        x = x->Next;
        delete y;
    }

    First = NULL;
    Last = NULL;
    Root = NIL;
    mCount = 0;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::GetFirst() const
{
    return First;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::GetLast() const
{
    return Last;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
int             RBTree<K,D,C>::GetCount()
{
    return mCount;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::GetNextUnique(RBTreeNode<K,D,C> * node)
{
    assert(node);

    if (node->Color == scDUP) node = Find(node->Key);
    if (node) node = TreeNext(node);

    return node;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::GetPrevUnique(RBTreeNode<K,D,C> * node)
{
    assert(node);

    if (node->Color == scDUP) node = Find(node->Key);
    if (node) node = TreePrev(node);

    return node;
}


//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::RotateRight(RBTreeNode<K,D,C> * x)
{
    RBTreeNode<K,D,C> * y = x->Left;

    x->Left = y->Right;
    if (y->Right != NIL) y->Right->Parent = x;

    if (y != NIL) y->Parent = x->Parent;
    if (x->Parent)
    {
        if (x == x->Parent->Right)  x->Parent->Right = y;
        else                        x->Parent->Left = y;
    }
    else
    {
        Root = y;
    }

    y->Right = x;
    if (x != NIL) x->Parent = y;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::RotateLeft(RBTreeNode<K,D,C> * x)
{
    RBTreeNode<K,D,C> * y = x->Right;

    x->Right = y->Left;
    if (y->Left != NIL) y->Left->Parent = x;

    if (y != NIL) y->Parent = x->Parent;
    if (x->Parent)
    {
        if (x == x->Parent->Left)   x->Parent->Left = y;
        else                        x->Parent->Right = y;
    }
    else
    {
        Root = y;
    }

    y->Left = x;
    if (x != NIL) x->Parent = y;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::FixInsert(RBTreeNode<K,D,C> * x)
{
    while (x != Root && x->Parent->Color == scRED)
    {
        if (x->Parent == x->Parent->Parent->Left)
        {
            RBTreeNode<K,D,C> *  y = x->Parent->Parent->Right;
            if (y->Color == scRED)
            {
                x->Parent->Color = scBLACK;
                y->Color = scBLACK;
                x->Parent->Parent->Color = scRED;
                x = x->Parent->Parent;
            }
            else
            {
                if (x == x->Parent->Right)
                {
                    x = x->Parent;
                    RotateLeft(x);
                }

                x->Parent->Color = scBLACK;
                x->Parent->Parent->Color = scRED;
                RotateRight(x->Parent->Parent);
            }
        }
        else
        {
            RBTreeNode<K,D,C> *  y = x->Parent->Parent->Left;
            if (y->Color == scRED)
            {
                x->Parent->Color = scBLACK;
                y->Color = scBLACK;
                x->Parent->Parent->Color = scRED;
                x = x->Parent->Parent;
            }
            else
            {
                if (x == x->Parent->Left)
                {
                    x = x->Parent;
                    RotateRight(x);
                }
                x->Parent->Color = scBLACK;
                x->Parent->Parent->Color = scRED;
                RotateLeft(x->Parent->Parent);
            }
        }
    }

    Root->Color = scBLACK;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::DoInsert(RBTreeNode<K,D,C> * x)
{
    RBTreeNode<K,D,C> *  y;
    RBTreeNode<K,D,C> *  p;
    int c=0;

    y = Root;
    p = NULL;

    while (y != NIL)
    {
        c = ComparFn(x->Key, y->Key);
        if (c == 0)
        {
            x->Color = scDUP;
            y = TreeNext(y);
            if (y)  DLLInsertBefore(y, x);
            else    DLLAddAfter(Last, x);

            return;
        }

        p = y;
        y = (c < 0) ? y->Left : y->Right;
    }

    if (p)
    {
        x->Parent = p;

        if (c < 0)
        {
            DLLInsertBefore(p, x);

            p->Left = x;
        }
        else
        {
            y = TreeNext(p);
            if (y)  DLLInsertBefore(y, x);
            else    DLLAddAfter(Last, x);

            p->Right = x;
        }
    }
    else
    {
        Root = x;
        x->Color = scBLACK;
        First = x;
        Last = x;
    }

    FixInsert(x);
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::FixDelete(RBTreeNode<K,D,C> * x)
{
    while (x != Root && x->Color == scBLACK)
    {
        if (x == x->Parent->Left)
        {
            RBTreeNode<K,D,C> * w = x->Parent->Right;
            if (w->Color == scRED)
            {
                w->Color = scBLACK;
                x->Parent->Color = scRED;
                RotateLeft (x->Parent);
                w = x->Parent->Right;
            }
            if (w->Left->Color == scBLACK && w->Right->Color == scBLACK)
            {
                w->Color = scRED;
                x = x->Parent;
            }
            else
            {
                if (w->Right->Color == scBLACK)
                {
                    w->Left->Color = scBLACK;
                    w->Color = scRED;
                    RotateRight(w);
                    w = x->Parent->Right;
                }
                w->Color = x->Parent->Color;
                x->Parent->Color = scBLACK;
                w->Right->Color = scBLACK;
                RotateLeft (x->Parent);
                x = Root;
            }
        }
        else
        {
            RBTreeNode<K,D,C> * w = x->Parent->Left;
            if (w->Color == scRED) {
                w->Color = scBLACK;
                x->Parent->Color = scRED;
                RotateRight (x->Parent);
                w = x->Parent->Left;
            }
            if (w->Right->Color == scBLACK && w->Left->Color == scBLACK)
            {
                w->Color = scRED;
                x = x->Parent;
            }
            else
            {
                if (w->Left->Color == scBLACK)
                {
                    w->Right->Color = scBLACK;
                    w->Color = scRED;
                    RotateLeft (w);
                    w = x->Parent->Left;
                }
                w->Color = x->Parent->Color;
                x->Parent->Color = scBLACK;
                w->Left->Color = scBLACK;
                RotateRight (x->Parent);
                x = Root;
            }
        }
    }
    x->Color = scBLACK;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::DoDelete(RBTreeNode<K,D,C> * x)
{
    if (x->Color == scDUP)
    {
        DLLDelete(x);
        return;
    }

    if (x->Next && x->Next->Color == scDUP)
    {
        ReplaceTreeNode(x, x->Next);
        DLLDelete(x);
        return;
    }

    DLLDelete(x);

    RBTreeNode<K,D,C> *  y;
    RBTreeNode<K,D,C> *  z;
    bool            balance;

    if (x->Left == NIL || x->Right == NIL)  y = x;
    else                                    y = FindMin(x->Right);

    if (y->Left != NIL) z = y->Left;
    else                z = y->Right;

    z->Parent = y->Parent;
    if (y->Parent)
    {
        if (y == y->Parent->Left)   y->Parent->Left = z;
        else                        y->Parent->Right = z;
    }
    else
    {
        Root = z;
    }

    balance = (y->Color == scBLACK);

    if (y != x) ReplaceTreeNode(x, y);

    if (balance) FixDelete(z);
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::DLLDelete(RBTreeNode<K,D,C> * x)
{
    if (x->Prev)    x->Prev->Next = x->Next;
    else            First = x->Next;
    if (x->Next)    x->Next->Prev = x->Prev;
    else            Last = x->Prev;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::DLLInsertBefore(RBTreeNode<K,D,C> * x, RBTreeNode<K,D,C> * y)
{
    y->Prev = x->Prev;
    if (y->Prev)    y->Prev->Next = y;
    else            First = y;
    y->Next = x;
    x->Prev = y;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::DLLAddAfter(RBTreeNode<K,D,C> * x, RBTreeNode<K,D,C> * y)
{
    y->Next = x->Next;
    if (y->Next)    y->Next->Prev = y;
    else            Last = y;
    y->Prev = x;
    x->Next = y;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::InitNode(RBTreeNode<K,D,C> * x)
{
    x->Next = NULL;
    x->Prev = NULL;
    x->Parent = NULL;
    x->Left = NIL;
    x->Right = NIL;
    x->Color = scRED;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::TreeNext(RBTreeNode<K,D,C> * x)
{
    if (x->Right != NIL) return FindMin(x->Right);

    while (x->Parent)
    {
        if (x == x->Parent->Left) break;

        x = x->Parent;
    }

    return x->Parent;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::TreePrev(RBTreeNode<K,D,C> * x)
{
    if (x->Left != NIL) return FindMax(x->Left);

    while (x->Parent)
    {
        if (x == x->Parent->Right) break;

        x = x->Parent;
    }

    return x->Parent;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::FindMin(RBTreeNode<K,D,C> * x)
{
    while (x->Left != NIL) x = x->Left;
    return x;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
RBTreeNode<K,D,C> *  RBTree<K,D,C>::FindMax(RBTreeNode<K,D,C> * x)
{
    while (x->Right != NIL) x = x->Right;
    return x;
}

//---------------------------------------------------------------------------
template <class K, class D, typename C>
inline
void            RBTree<K,D,C>::ReplaceTreeNode(RBTreeNode<K,D,C> * x, RBTreeNode<K,D,C> * y)
{
    y->Color = x->Color;
    y->Left = x->Left;
    y->Right = x->Right;
    y->Parent = x->Parent;

    if (y->Parent)
    {
        if (y->Parent->Left == x)   y->Parent->Left = y;
        else                        y->Parent->Right = y;
    }
    else
    {
        Root = y;
    }

    y->Left->Parent = y;
    y->Right->Parent = y;
}

}

//---------------------------------------------------------------------------
#endif

