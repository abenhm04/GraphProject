package graphexplorer;

/**
 * A test bed for a weighted digraph abstract data type implementation
 * and implementations of elementary classical graph algorithms that use the
 * ADT
 * @see GraphAPI.java, Graph.java, City.java
 * @author Duncan, Adam Ben Hmida
 * <pre>
 * usage: GraphExplorer <graphFileName>
 * <graphFileName> - a text file containing the description of the graph
 * in DIMACS file format
 * Date: 11-20-2023
 * course: csc 3102
 * programming project 3
 * Instructor: Dr. Duncan
 * </pre>
 */

import java.io.*;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.Comparator;
import java.util.PriorityQueue; 
import java.util.concurrent.atomic.AtomicBoolean;

public class GraphExplorer
{
   public static final Double INFINITY = Double.POSITIVE_INFINITY;
   public static final Integer NIL = -1;

   public static void main(String []args) throws GraphException
   {
      if (args.length != 1)
      {
         System.out.println("Usage: GraphExplorer <filename>");
         System.exit(1);
      }
      City city1, city2;
      Scanner console = new Scanner(System.in);
      int menuReturnValue, i,j;
      //Define a lambda function f:City -> PrintStream that displays
      //the key followed by the label (name of the city) in two-column format
      //and advances the cursor to the next line. This function will be used
      //as an argument to bfsTraverse and dfsTraverse and to display the 
      //topological listing of the vertices.
      //Function<City,PrintStream> f = ...;	  
      Function<City, Void> f = city -> {
    	    System.out.printf("%-5d %-20s\n", city.getKey(), city.getLabel());
    	    return null; // Return type is Void, so return null
    	};
    	
      Graph<City> g = readGraph(args[0]);      
      long s = g.size();
      menuReturnValue = -1;
      while (menuReturnValue != 0)
      {
         menuReturnValue = menu();
         switch(menuReturnValue)
         {
            case 1: //{-1,0,1} Incidence Matrix of g
               System.out.println();
               int[][] iMat = genIncMatUDG(g);
               if (iMat == null)
               {
                   System.out.println("The graph in "+args[0]+" has no edges.");
               }
               else
               {
                   System.out.println("|V| x |E| {-1,0,1}-Incidence Matrix of g in "+args[0]);
                   System.out.println("--------------------------------------------------------");
                   for (i = 0; i < g.size(); i++)
                   {
                       for (j = 0; j < g.countEdges(); j++)
                       {
                           System.out.printf("%3d",iMat[i][j]);
                       }
                       System.out.println();
                   }
               }
               System.out.println();
               break;
            case 2: //post-order depth-first-search traversal of g'
               System.out.println();
               System.out.println("Breadth-First-Search Traversal of the Graph In "+args[0]);
               System.out.println("===============================================================");
               //invoke the bfsTraverse function
               // Output should be aligned in two-column format as illustrated below:
               // 1     Charlottetown
               // 4     Halifax
               // 2     Edmonton     
               g.bfsTraverse(f); 

               System.out.println("===============================================================");
               System.out.println();
               System.out.println("Depth-First-Search Traversal of the Graph In "+args[0]);
               System.out.println("===============================================================");
               //invoke the dfsTraverse function
               // Output should be aligned in two-column format as illustrated below:
               // 1     Charlottetown
               // 4     Halifax
               // 2     Edmonton     
               g.dfsTraverse(f);

               System.out.println("===============================================================");
               System.out.println();
               System.out.println();
               break;
            case 3: //Check Strong Connectivity of G
                 System.out.println();
                 System.out.printf("The digraph %s is ",args[0]);
                 if (isStronglyConnected(g))
                 {
                     System.out.printf("strongly connected.%n");
                 }
                 else
                 {
                     System.out.printf("not strongly connected.%n");
                 }
                 System.out.println();
                 System.out.println();
                 break;                            
            case 4://out-degree topological sort
                 System.out.println();
                 int[] top = new int[(int)g.size()];
                 if (topSortOutDeg(g,top))
                 {
                    System.out.println("Topological Sorting of The Graph In "+args[0]);
                    System.out.println("===============================================================");           

                    for (int vertexIndex = 0; vertexIndex < top.length; vertexIndex++) {
                        City city = g.retrieveVertex(new City(top[vertexIndex]));
                        System.out.printf("%-5d %-20s\n", city.getKey(), city.getLabel().trim());
                    }
                    System.out.println("===============================================================");
                 }
                 else
                    System.out.println("No topological ordering possible. The digraph in "+args[0]+" contains a directed cycle.");
                 System.out.printf("%n%n");
                 top = null;
                 break;
            case 5: //KruskalMST;
               System.out.printf("%nThe Minimum Spanning Tree of the graph in "+args[0]+".%n%n");
               int edgesInMST = 0;
               int[] mstK = new int[(int)g.size()];
               double totalWeight = kruskalMST(g,mstK);
               String cityNameA,cityNameB;
               for (i=1; i<=g.size(); i++)
               {
                   if (mstK[i-1] == -1)
                       cityNameA = "NONE";
                   else
                   {
                       edgesInMST++;
                       cityNameA = g.retrieveVertex(new City(mstK[i-1])).getLabel().trim();
                   }
                   cityNameB = g.retrieveVertex(new City(i)).getLabel().trim();
                   System.out.printf("%d-%s parent[%d] <- %d (%s)%n",i,cityNameB,i,mstK[i-1],cityNameA);
               }
               System.out.printf("%nThe weight of the minimum spanning tree/forest is %.2f miles.%n",totalWeight);    
			   System.out.printf("Spanning Tree Edge Density is %.2f%%.%n%n",100.0*edgesInMST/g.countEdges());
               break;
            default: ;
         } //end switch
      }//end while
   }//end main

   /**
    * This method reads a text file formatted as described in the project description.
    * @param filename the name of the DIMACS formatted graph file.
    * @return an instance of a graph.
    */
   private static Graph<City> readGraph(String filename)
   {
      try
      {
         Graph<City> newGraph = new Graph();
         try (FileReader reader = new FileReader(filename)) 
         {
            char temp;
            City city1, city2, aCity;
            String tmp;
            int k, m, v1, v2, j, size=0, nEdges=0;
            Integer key, v1Key,v2Key;
            Double weight;
            Scanner in = new Scanner(reader);
            while (in.hasNext())
            {
                tmp = in.next();
                temp = tmp.charAt(0);
                if (temp == 'p')
                {
                    size = in.nextInt();
                    nEdges = in.nextInt();
                }
                else if (temp == 'c')
                {
                    in.nextLine();
                }
                else if (temp == 'n')
                {
                    key = in.nextInt();
                    tmp = in.nextLine();
                    aCity = new City(key,tmp);
                    newGraph.insertVertex(aCity); 
                }
                else if (temp == 'e')
                {
                    v1Key = in.nextInt();
                    v2Key = in.nextInt();
                    weight = in.nextDouble();
                    city1 = new City(v1Key);
                    city2 = new City(v2Key);
                    newGraph.insertEdge(city1,city2,weight); 
                }
            }
         }
         return newGraph;
      }
      catch(IOException exception)
      {
            System.out.println("Error processing file: " + exception);
      }
      return null;
   } 

   /**
    * Display the menu interface for the application.
    * @return the menu option selected.
    */  
   private static int menu()
   {
      Scanner console = new Scanner(System.in);
      String option;
      do
      {
         System.out.println("  BASIC WEIGHTED GRAPH APPLICATION   ");
         System.out.println("=======================================================");
         System.out.println("[1] Generate {-1,0,1} Incidence Matrix of G");
         System.out.println("[2] BFS/DFS Traversal of G");
         System.out.println("[3] Check Strong Connectedness of G");
         System.out.println("[4] Topological Ordering of V(G)");
         System.out.println("[5] Kruskal's Minimum Spanning Tree/Forest in G");
         System.out.println("[0] Quit");
         System.out.println("=====================================");
         System.out.printf("Select an option: ");         
         option = console.nextLine().trim();
         try
         {
             int choice = Integer.parseInt(option);
             if (choice < 0 || choice > 6)
             {
                System.out.println("Invalid option...Try again");
                System.out.println();
             }
             else
                return choice;
         }
         catch(NumberFormatException e)
         {
            System.out.println("Invalid option...Try again");
         }                           
      }while(true);
   }
   
   /**
    * Generates a {-1,0,1} vertices x edges incidence matrice for a
    * digraph; the columns representing the edges of the graph are arranged 
    * in lexicographical order, source=1 -> destination=-1 (source,destination).
    * @param g a digraph
    * @return a |V| x |E| 2-D array representing the {-1,0,1} incidence matrix
    * of g such that m[i-1,e-1] = 1 and m[j-1,e-1] = -1 if there is an outgoing
    * edge from vertex i to vertex j representing the eth edge (i->j) in g when the
    * edges are arranged in lexicographical order (source,destination); all otherwise
    * entries in a column is 0 when there isn't an incident edge. The null reference is
    * returned when g contains no edges.
    */
   private static int[][] genIncMatUDG(Graph<City> g)  throws GraphException
   {
          if (g == null || g.countEdges() == 0) {
               return null;
           }

           int numVertices = (int) g.size();
           int numEdges = (int) g.countEdges();

           // Initialize the incidence matrix with zeros
           int[][] incMat = new int[numVertices][numEdges];

           // Edge counter
           int edgeCounter = 0;

           // Iterate through all pairs of vertices to populate the incidence matrix
           for (int i = 1; i <= numVertices; i++) {
               for (int j = 1; j <= numVertices; j++) {
                   if (g.isEdge(new City(i), new City(j))) {
                       boolean Undirected = g.isEdge(new City(j), new City(i));

                       incMat[i - 1][edgeCounter] = 1;
                       incMat[j - 1][edgeCounter] = Undirected ? 1 : -1;

                       edgeCounter++;
                   }
               }
           }

           return incMat;
       }

		

   
   /**
    * Determines whether the specified digraph is strongly connected
    * @param g a digraph
    * @return true if the digraph is strongly connected; otherwise, false
    */
   private static boolean isStronglyConnected(Graph<City> g)  throws GraphException
   {
	   //implement this method
	   if (g == null || g.size() == 0) {
	        return false;
	    }

	    int numVertices = (int) g.size();

	    // Check reachability from each vertex
	    for (int i = 1; i <= numVertices; i++) {
	        Set<Integer> visited = new HashSet<>();
	        Stack<Integer> stack = new Stack<>();

	        // Initialize DFS from vertex i
	        stack.push(i);
	        visited.add(i);

	        while (!stack.isEmpty()) {
	            int current = stack.pop();

	            // Check all possible vertices to see if they're adjacent (connected)
	            for (int adjKey = 1; adjKey <= numVertices; adjKey++) {
	                if (g.isEdge(new City(current), new City(adjKey)) && !visited.contains(adjKey)) {
	                    stack.push(adjKey);
	                    visited.add(adjKey);
	                }
	            }
	        }

	        if (visited.size() != numVertices) {
	            // If not all vertices are reachable from this vertex, the graph is not strongly connected
	            return false;
	        }
	    }

	    return true;
   }

   /**
    * Merges two trees in a forest where (i,j) is a bridge between
    * trees.
    * @param i a vertex in one tree (source)
    * @param j a vertex in another tree (destination)
    * @param parent the parent-array representation of a forest
    * where the parent of a root is -1
    */
   public static void trUnion(int i, int j, int[] parent) {
		Stack<Integer> stackJ = new Stack<>();
		if(find(parent,i)!= find(parent,j)) {
			int rootJ = j;
			while (parent[rootJ-1] != -1) {
				stackJ.push(rootJ);
				rootJ = parent[rootJ-1];
			}
			stackJ.push(rootJ);
		
		
		while (stackJ.size()>1) {
			int tempValue = stackJ.pop();
			parent[tempValue - 1] = stackJ.peek();
		}
		
		int tempValue = stackJ.pop();
		parent[tempValue-1] = i;
		
		} 
	}
   /**
    * Finds the root vertex of the tree in which the specified is
    * @param parent the parent implementation of disjointed subtrees of a graph
    * @param v a vertex
    * @return the root of the subtree containing v or v if v is a root.
    */  
   private static int find(int[] parent, int v) {
		for (int i = 0; i< parent.length; i++) {
			if (parent[v-1] == -1) {
				return v;
			}
			else if (i+1 == v) {
				return find(parent, parent[i]);
			}
		}
		return -1;
	}

   /**
    * This method generates a minimum spanning tree using Kruskal's 
    * algorithm. If no such MST exists, then it generates a minimum spanning forest.
    * @param g a weighted directed graph
    * @param parent the parent implementation of the minimum spanning tree/forest
    * @return the weight of such a tree or forest.
    * @throws GraphException when this graph is empty
    * <pre>
    * {@code
    * If a minimum spanning tree cannot be generated,
    * the parent implementation of a minimum spanning tree or forest is
    * determined. This implementation is based on the union-find stragey.
    * }
    * </pre>
    */    
   private static double kruskalMST(Graph<City> g, int[] parent) throws GraphException {
		/**
		 * An auxiliary data structure to store the information about an edge for the
		 * MST algorithm
		 * 
		 */
		class EdgeType {
			public double weight;
			public int source;
			public int destination;
			public boolean chosen;
		}
		/* An EdgeType comparator */
		Comparator<EdgeType> cmp = (t1, t2) -> {
			if (t1.weight < t2.weight)
				return -1;
			if (t1.weight > t2.weight)
				return 1;
			if (t1.source < t2.source)
				return -1;
			if (t1.source > t2.source)
				return 1;
			if (t1.destination < t2.destination)
				return -1;
			if (t1.destination > t2.destination)
				return 1;
			return 0;
		};
		// implement this method
		/*
		 * Defining an instance of the PriorityQueue class that uses the comparator
		 * above and complete the implementation of the algorithm
		 */

		if (g == null || g.size() == 0) {
			return 0.0;
		}

		int numberVertices = (int) g.size();
		PriorityQueue<EdgeType> pq = new PriorityQueue<>(cmp);
		for (int i = 1; i <= numberVertices; i++) {
			parent[i - 1] = -1;
		}
		for (int i = 1; i <= numberVertices; i++) {
		    for (int j = 1; j <= numberVertices; j++) { 
		        if (g.isEdge(new City(i), new City(j))) { 
		            EdgeType edge = new EdgeType();
		            edge.source = i;
		            edge.destination = j;
		            edge.weight = g.retrieveEdge(new City(i), new City(j));
		            pq.add(edge);
		        }
		    }
		}
		
		double totalWeight = 0.0;
		int count = 0;
		int v = (int)g.size()-1;
		while (!pq.isEmpty()&& count < v) {
			EdgeType edge = pq.poll();
			int root1 = find(parent, edge.source);
			int root2 = find(parent, edge.destination);
			if (root1 != root2) {
				trUnion(edge.source, edge.destination, parent);
				totalWeight += edge.weight;
				count++;
			}
		}
		return totalWeight;
	}

            
            
         
   /**
    * Generates a topological labeling of the specified digraph, in reverse order,
    * using the decrease-and-conquer algorithm that successively selects and r 
    * removes a vertex of out-degree 0 until all the vertices are selected. 
    * The graph is explored in lexicographical order when adding a new vertex to the  
    * topological ordering and the graph is not modified. Updates of the degrees 
    * and vertices that are selected are tracked using auxiliary data structures. 
    * @param g a digraph
    * @param linearOrder the topological ordering of the vertices 
    * @return true if a topological ordering of the vertices of the specified digraph 
    * exists; otherwise, false. 
    */   
       private static boolean topSortOutDeg(Graph<City> g, int linearOrder[]) throws GraphException
       {
           if (g == null || g.isEmpty()) {
               return false;
           }

           int numVertices = (int)g.size();

           Map<Integer, Integer> outDegMap = new HashMap<>();
           for (int i = 1; i <= numVertices; i++) {
               City city = new City(i);
               outDegMap.put(i, (int) g.outDegree(city));

           }

           int exitVertex = -1;
           int tmpDeg;
           int index = numVertices-1;

           for(int i = 1; i<= numVertices; i++) {
               tmpDeg = outDegMap.get(i);
               if(tmpDeg == 0) {
                   exitVertex = i;
                   linearOrder[index] = exitVertex;
                   index--;
                   for (int j = 1; j <= numVertices; j++) {
                        City adjacentCity = g.retrieveVertex(new City(j));
                        if (g.isEdge(adjacentCity, new City(exitVertex))) {
                            outDegMap.put(j, outDegMap.get(j) - 1);
                        }
                    }
                   outDegMap.put(exitVertex, -1);
                   i = 0;

               }
           }
           if(index>0) {
               return false;
           }

           for(int i=0; i<linearOrder.length; i++) {
               City tempcity = g.retrieveVertex(new City(linearOrder[i]));
               System.out.println(linearOrder[i] + "\t" + tempcity.getLabel().trim());
           }

           if(exitVertex == -1) {
               return false;
           }

           return true;
         }
}
