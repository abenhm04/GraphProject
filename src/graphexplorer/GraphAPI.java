package graphexplorer;
//package graphexplorer;
import java.util.function.Function;

/**
 * Describe the fundamental operations of a weighted directed graph and defines
 * a class for reporting exceptions on operations of the abstract data type.
 * @author William Duncan
 * <pre>
 * Date: 99-99-9999
 * Course: csc 3102
 * Programming Project # 3 
 * Instructor: Dr. Duncan
 * Note: DO NOT MODIFY THIS FILE 
 * </pre>
 */

/**
 * For reporting Graph-related exception
 * @author Duncan
 */
class GraphException extends Exception 
{
    /**
     * Creates a new instance of <code>GraphException</code> without detail
     * message.
     */
    public GraphException() 
    {
    }

    /**
     * Constructs an instance of <code>GraphException</code> with the specified
     * detail message.
     * @param msg the detail message.
     */
    public GraphException(String msg) 
    {
        super(msg);
    }
}

/**
 * Describes the fundamental operations of a weighted directed graph.
 * @param <E> the data type
 * @author William Duncan
 * @param <E> The element type of the graph ADT
 */
public interface GraphAPI<E extends Comparable>
{
   /**
    * This method inserts a new vertex whose data item is data in the
    * weighted graph. If the key of the data already exists the vertex
    * is updated. No two vertices may have the same key.
    * @param data - data stored in a vertex.
    */
    void insertVertex(E data);

   /**
    * This method deletes a vertex whose data item is data from the
    * weighted graph. If the key of the data does not exist or the 
    * in-degree or out degree is positive, the weighted graph remains 
    * unaltered. It the data item exists and the in-degree
    * and out-degree are both 0, the vertex is removed from the graph.
    * @param key - data stored in a vertex.
    */
    void deleteVertex(E key);

   /**
    * This method inserts a weighted directed edge between two vertices.
    * If either key does not exist, the graph remains unaltered. If there
    * is already a directed edge between the vertices, its weight is
    * updated. If both keys exist and there isn't an edge between the
    * vertices, an edge is inserted and the in-degree and out-degree of 
    * both vertices are updated.
    * @param fromKey - data of the originating vertex.
    * @param toKey - data of the destination vertex.
    * @param weight - weight of the edge between the from and to
    * vertices.
    */
    void insertEdge(E fromKey, E toKey, Double weight);

   /**
    * This method removes a weighted directed edge between two vertices.
    * If either key does not exist, the method returns and the graph
    * remains unaltered. If the edge exists, it is removed from the vertex
    * with fromKey to the vertex with toKey. The in-degree and out-degree
    * of both vertices are updated. 
    * @param fromKey - search key of the originating vertex.
    * @param toKey - search key of the destination vertex.
    */
    void deleteEdge(E fromKey, E toKey);

   /**
    * This method returns the weight of the directed edge between
    * the vertices with fromKey and toKey if the edge exists. If
    * the directed edge does not exist, an exception is generated.
    * @param fromKey - search key of the originating vertex.
    * @param toKey - search key of the destination vertex.
    * @return the weight on the edge
    * @throws GraphException
    */
    double retrieveEdge(E fromKey, E toKey) throws GraphException;

   /** 
    * This method returns the item stored in the vertex whose key
    * is key. If the key does not exist, an exception is generated.
    * @param key - search key of the vertex.
    * @return 
    * @throws GraphException
    */
    E retrieveVertex(E key) throws GraphException;

   /**
    * This method applies the visit function to the vertices of
    * the graph in breadth-first-search order.
    * @param func - the visit function.
    */
    void bfsTraverse(Function func);

   /**
    * This method applies the visit function to the vertices of
    * the graph in postorder depth-first-search order.
    * @param func - the visit function.
    */
    void dfsTraverse(Function func);

   /**
    * Determine whether the graph is empty.
    * @return this function returns true if the graph is empty;
    * otherwise, it returns false if the graph contains at least one node.
    */
   boolean isEmpty();

   /**
    * Returns the order of graph.
    * @return the number of vertices in the graph
    */
   long size();

   /**
    * Determine whether an item is in the graph.
    * @param key search key of the vertex.
    * @return true on success; false on failure.
    */
   boolean isVertex(E key);

   /**
    * Determines whether there is a directed edge between two vertices.
    * @return true on success or false on failure.
    * @param fromKey - search key of the originating vertex.
    * @param toKey - search key of the destination vertex.
    */
   boolean isEdge(E fromKey, E toKey);
   
   /**
    * Determines the number of edges in the graph.
    * @return the number of edges.
    */
   long countEdges();

   /**
    * Determines the number of out-directed edges from the vertex
    * with the key.
    * @param key
    * @return out-degree.
    * @throws GraphException when the vertex with the specified
    * key does not exist
    */
   long outDegree(E key) throws GraphException;

   /**
    * Determines the number of in-directed edges from the vertex
    * with the key.
    * @param key
    * @return in-degree.
    * @throws GraphException when the vertex with the specified key
    * does not exist
    */
   long inDegree(E key) throws GraphException;
}

