package graphexplorer;
//package graphexplorer;

import java.util.Comparator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;
import java.util.function.Function;

/** 
 * Implements an interface for a weighted digraph
 * @param <E> the data type of the graph; E must implement that standard Java
 * parametric Comparable interface
 * @author Duncan, Adam Ben Hmida
 * @see GraphAPI
 * <pre>
 * Date: 11-21-2023
 * CSC 3102 Programming Project # 3
 * Instructor: Dr. Duncan 
 *
 * DO NOT REMOVE THIS NOTICE (GNU GPL V2):
 * Contact Information: duncanw@lsu.edu
 * Copyright (c) 2023 William E. Duncan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>
 * </pre>
 */
public class Graph<E extends Comparable<E>> implements GraphAPI<E>
{
    /*
     * number of vertices (size of this graph)
     */
    private long order;
    /**
     * pointer to the list of vertices
     */
    private Vertex first;
    
   /**
    * A comparator lambda function that compares two nodes (elements)
    * of this graph; cmp.compare(x,y) gives 1. negative when x less than y
    * 2. positive when x greater than y 3. 0 when x equal y
    */   
   private Comparator<? super E> cmp;          
   
   /**
    * A vertex of a graph stores a data item and references
    * to its edge list and the succeeding vertex. The data
    * object extends the comparable interface.
    */      
   private class Vertex
   {
       /**
        * pointer to the next vertex
        */
       public Vertex pNextVertex;
       /**
        * the data item
        */
       public E data;
       /**
        * in-degree
        */
       public long inDeg;
       /**
        * out-degree
        */
       public long outDeg;
       /**
        * pointer to the edge list
        */
      public Edge pEdge;
      /**
       * Field for tracking vertex accesses
       */
      public long processed;
   }
   
   /** 
    * An edge of a graph contains a reference to the destination
    * vertex, a reference to the succeeding edge in the edge list and
    * the weight of the directed edge.
    */   
   private class Edge
   {
       /**
        * pointer to the destination vertex
        */
       public Vertex destination;
       /**
        * weight on this edge
        */
       public Double weight;
       /**
        * pointer to the next edge
        */
       public Edge pNextEdge;
   }   
   /**
    * Constructs an empty weighted directed graph
   */
   public Graph()
   {
      first = null;
      order = 0;
      cmp = (x,y) -> x.compareTo(y);         
   }
   
   /**
    * A parameterized constructor that uses an externally defined comparator    
    * @param fn - a trichotomous integer value comparator function   
    */
   public Graph(Comparator<? super E> fn)
   {
      first = null;
      order = 0;
      cmp = fn;        
   }      

   @Override
   public void insertVertex(E obj)
   {
      Vertex locPtr = first;
      Vertex predPtr = null;
      while (locPtr != null && cmp.compare(obj,locPtr.data) > 0)
      {
         predPtr = locPtr;
         locPtr = locPtr.pNextVertex;
      }
      /*key already exist. */
      if (locPtr != null && cmp.compare(obj,locPtr.data)==0)
      {
         locPtr.data = obj;
         return;
      }
      Vertex newPtr = new Vertex();
      newPtr.pNextVertex = null;
      newPtr.data  = obj;
      newPtr.inDeg = 0;    
      newPtr.outDeg = 0;    
      newPtr.processed = 0;
      newPtr.pEdge = null;      
      /* insert before first vertex */
      if (predPtr == null)
         first = newPtr;
      else
         predPtr.pNextVertex = newPtr;
      newPtr.pNextVertex = locPtr;   
      order++;
   }

   @Override
   public void deleteVertex(E key)
   {
      if (isEmpty())
         return;
      Vertex predPtr = null;
      Vertex walkPtr = first;
      while (walkPtr != null && cmp.compare(key,walkPtr.data) > 0)
      {
         predPtr = walkPtr;
         walkPtr = walkPtr.pNextVertex;
      }
      if (walkPtr == null || cmp.compare(key,walkPtr.data) != 0)
         return;
      /* vertex found. Test degree */
      if ((walkPtr.inDeg > 0) || (walkPtr.outDeg > 0))
         return;
      /* OK to delete */
      if (predPtr == null)
         first = walkPtr.pNextVertex;
      else
         predPtr.pNextVertex = walkPtr.pNextVertex;
      order--;
   }

   @Override
   public void insertEdge(E fromKey, E toKey, Double weight)
   {      
      if (isEmpty())
         return;
      Edge pred;   
      /* check whether originating vertex exists */
      Vertex tmpFrom = first;
      while (tmpFrom != null && cmp.compare(fromKey,tmpFrom.data) > 0)
         tmpFrom = tmpFrom.pNextVertex;
      if (tmpFrom == null || cmp.compare(fromKey,tmpFrom.data) != 0)
         return;
      /* locate destination vertex */
      Vertex tmpTo = first;
      while (tmpTo != null && cmp.compare(toKey,tmpTo.data) > 0)
         tmpTo = tmpTo.pNextVertex;
      if (tmpTo == null || cmp.compare(toKey,tmpTo.data) != 0)
         return;
      /*check if edge already exists */
      Edge tmpEdge = tmpFrom.pEdge;
      while (tmpEdge != null && tmpEdge.destination != tmpTo)
         tmpEdge = tmpEdge.pNextEdge;
      if (tmpEdge != null && tmpEdge.destination != null)
         return;
      tmpFrom.outDeg++;
      tmpTo.inDeg++;
      Edge newEdge = new Edge();
      newEdge.destination = tmpTo;
      newEdge.weight = weight;
      newEdge.pNextEdge = null;
      if (tmpFrom.pEdge == null)
      {
         tmpFrom.pEdge = newEdge;
         return;
      }
      pred = null;
      tmpEdge = tmpFrom.pEdge;
      while (tmpEdge != null && tmpEdge.destination != tmpTo
	  &&cmp.compare(newEdge.destination.data,tmpEdge.destination.data)>0)
      {
         pred = tmpEdge;
         tmpEdge = tmpEdge.pNextEdge;
      }
      if (pred == null)
         tmpFrom.pEdge = newEdge;
      else
         pred.pNextEdge = newEdge;
      newEdge.pNextEdge = tmpEdge;
   }         

   @Override
   public void deleteEdge(E fromKey, E toKey)
   {
      /* find source vertex */
      Vertex tmpFrom = first;
      while (tmpFrom != null && cmp.compare(fromKey,tmpFrom.data)>0)
         tmpFrom = tmpFrom.pNextVertex;
      if (tmpFrom == null || cmp.compare(fromKey,tmpFrom.data) != 0)
         return;
      /* locate destination vertex */
      Vertex tmpTo = first;
      while (tmpTo != null && cmp.compare(toKey,tmpTo.data)>0)
         tmpTo = tmpTo.pNextVertex;
      if (tmpTo == null || cmp.compare(toKey,tmpTo.data) != 0)
         return;
      /*check if edge does not exist */
      Edge tmpEdge = tmpFrom.pEdge;
      Edge pred = null;
      while (tmpEdge != null && tmpEdge.destination != tmpTo)
      {
         pred = tmpEdge;
         tmpEdge = tmpEdge.pNextEdge;
      }
      /* if edge does not exist */
      if (tmpEdge == null)
         return;
      /* update degrees */
      if (pred != null)
          pred.pNextEdge = tmpEdge.pNextEdge;
      tmpFrom.outDeg--;
      tmpTo.inDeg--;
   }

   @Override
   public double retrieveEdge(E fromKey, E toKey) throws GraphException
   {
      /* find source vertex */
      Vertex tmpFrom = first;
      while (tmpFrom != null && cmp.compare(fromKey,tmpFrom.data)>0)
         tmpFrom = tmpFrom.pNextVertex;
      if (tmpFrom == null || cmp.compare(fromKey,tmpFrom.data) != 0)
         throw new GraphException("Non-existent edge - retrieveEdge().");
      /* locate destination vertex */
      Vertex tmpTo = first;
      while (tmpTo != null && cmp.compare(toKey,tmpTo.data)>0)
         tmpTo = tmpTo.pNextVertex;
      if (tmpTo == null || cmp.compare(toKey,tmpTo.data) != 0)
         throw new GraphException("Non-existent edge - retrieveEdge().");
      /*check if edge does not exist */
      Edge tmpEdge = tmpFrom.pEdge;
      while (tmpEdge != null && tmpEdge.destination != tmpTo)
         tmpEdge = tmpEdge.pNextEdge;
      /* if edge does not exist */
      if (tmpEdge == null)
         throw new GraphException("Non-existent edge - retrieveEdge().");
      return tmpEdge.weight;
   }

   @Override
   public E retrieveVertex(E key) throws GraphException
   {
      if (isEmpty())
         throw new GraphException("Non-existent vertex - retrieveVertex().");
      Vertex tmp = first;
      while (tmp != null && cmp.compare(key,tmp.data) != 0)
         tmp = tmp.pNextVertex;
      if (tmp == null)
         throw new GraphException("Non-existent vertex - retrieveVertex().");
      return tmp.data;
   }

   @Override
   public void bfsTraverse(Function func)
   {
      if (isEmpty())
         return;
      Vertex walkPtr = first;
      while (walkPtr != null)
      {
         walkPtr.processed = 0;
         walkPtr = walkPtr.pNextVertex;
      }
      Queue<Vertex> queue = new LinkedList();
      Vertex toPtr;
      Edge edgeWalk;
      Vertex tmp;      
      walkPtr = first;
      while (walkPtr != null)
      {
         if (walkPtr.processed < 2)
         {
            if (walkPtr.processed < 1)
            {
               queue.add(walkPtr);
               walkPtr.processed = 1;
            }
         }
         while (!queue.isEmpty())
         {
            tmp = queue.remove();
            func.apply(tmp.data);
            tmp.processed = 2;
            edgeWalk = tmp.pEdge;
            while (edgeWalk != null)
            {
               toPtr = edgeWalk.destination;
               if (toPtr.processed == 0)
               {
                  toPtr.processed = 1;
                  queue.add(toPtr);
               }
               edgeWalk = edgeWalk.pNextEdge;
            }            
         }
         walkPtr = walkPtr.pNextVertex;      
      }
   }
   
   @Override
   public void dfsTraverse(Function func)
   {
      if (isEmpty())
         return;
      Vertex walkPtr = first;
      while(walkPtr != null)
      {
         walkPtr.processed = 0;
         walkPtr = walkPtr.pNextVertex;
      }
      Stack<Vertex> stack = new Stack();
      Vertex toPtr;      
      Edge edgeWalk;
      Vertex tmp;      
      walkPtr = first;
	  boolean hasUnexploredNbrs;
      while (walkPtr != null)
      {
         if (walkPtr.processed < 2)
         {
            if (walkPtr.processed < 1)
            {
               walkPtr.processed = 1;
               stack.push(walkPtr);
            }
            while (!stack.isEmpty())
            {
                tmp = stack.peek();
                edgeWalk = tmp.pEdge;
				hasUnexploredNbrs = false;
                while(edgeWalk != null)
                {
                   toPtr = edgeWalk.destination;
                   if (toPtr.processed == 0)
                   {
                      toPtr.processed = 1;
                      stack.push(toPtr);
                      edgeWalk = toPtr.pEdge;
					  hasUnexploredNbrs = true;
                   }
                   else
                      edgeWalk = edgeWalk.pNextEdge;
                }
				if (!hasUnexploredNbrs)
				{
				   tmp = stack.pop();
				   func.apply(tmp.data);
				   tmp.processed = 2;
				}
            }
         }
         walkPtr = walkPtr.pNextVertex;
      }            
   }

   @Override
   public boolean isEmpty()
   {
      return first == null;
   }

   @Override
   public long size()
   {
      return order;
   }

   @Override
   public boolean isVertex(E key)
   {
      if (isEmpty())
         return false;
      Vertex tmp = first;
      while (tmp != null && cmp.compare(key,tmp.data) != 0)
         tmp = tmp.pNextVertex;
      return tmp != null;      
   }

/***------------------  BEGIN AUGEMENTED METHOD ----------------***/
   @Override
   public boolean isEdge(E fromKey, E toKey)
   {
       if (isEmpty()) {
           return false;
       }

       try {
           double edgeWeight = retrieveEdge(fromKey, toKey);
           return true;
       } catch (GraphException e) {
           return false;
       }
   }

   @Override
   public long countEdges() {
       if (isEmpty()) {
           return 0;
       }

       long edgeCount = 0;
       Vertex vertex = first;

       while (vertex != null) {
           edgeCount += vertex.outDeg;
           vertex = vertex.pNextVertex;
       }

       return edgeCount;
   }

  
  
   
/***------------------  END AUGEMENTED METHOD ----------------***/

   @Override
   public long outDegree(E key) throws GraphException
   {
      if (isEmpty())
         throw new GraphException("Non-existent vertex - outDegree().");
      Vertex tmp = first;
      while (tmp != null && cmp.compare(key,tmp.data) != 0)
         tmp = tmp.pNextVertex;
      if (tmp == null)
         throw new GraphException("Non-existent vertex - outDegree().");
      return tmp.outDeg;
   }

   @Override
   public long inDegree(E key) throws GraphException
   {
      if (isEmpty())
         throw new GraphException("Non-existent vertex - inDegree().");
      Vertex tmp = first;
      while (tmp != null && cmp.compare(key,tmp.data) != 0)
         tmp = tmp.pNextVertex;
      if (tmp == null)
         throw new GraphException("Non-existent vertex - inDegree().");
      return tmp.inDeg;
   }    
}

