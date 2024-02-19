/**
 * CPA Mini projet 2023
 * 
 * ETUDIANTS STL :
 *  - Elhadj Alseiny DIALLO - 21314820
 *  - Tracy HONG - 21314944
 */

package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.Comparator;

public class DefaultTeam {
  static final int MIN_POINTS = 23 ;   // minimum number of points to keep after optimization       (after multiple tests, 23 is the best value)
  static final int BUDGET = 1664;       // maximum distance of the tree
  private Double maxOpti = Double.POSITIVE_INFINITY;
  private Point firstPoint;

  Comparator<Point> comp = new Comparator<Point>() {
    public int compare(Point a1, Point a2) {
      return (int) (a2.distance(firstPoint) - a1.distance(firstPoint));
    }
  };

  /**
   * Steiner algorithm.
   * This function is called in calculSteiner, calculSteinerBudget and optimiser.
   * @param points the list of points
   * @param edgeThreshold the maximum distance between two points to create an edge
   * @param hitPoints the list of hitPoints
   * @return ArrayList<Edge> the list of edges of the steiner tree
   */
  private ArrayList<Edge> steiner (ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints, int[][] paths) {
    ArrayList<Edge> edges = calculKruskal(hitPoints);
    ArrayList<Edge> edgesToKeep = new ArrayList<Edge>();

    //for each edge in kruskal, replace the edge by the shortest path in paths 
		for (Edge e : edges){
			ArrayList<Point> path = getPathPoints(points, paths, points.indexOf(e.p), points.indexOf(e.q)); 
			for (int l = 0; l < path.size()-1; l++){
				edgesToKeep.add(new Edge(path.get(l), path.get(l+1)));
			}
		}
    // System.out.println("edges: " + edgesToKeep.size()); // 32              //----------------------debug print
    // System.out.println("edgesToKeep: " + edgesToKeep.size()); // 73        //----------------------debug print

    //rebuild the list of points from the list of edges to keep 
    //and recalculate the minimum spanning tree to make sure the edges are valid (within the edgeThreshold)
    ArrayList<Point> newPoints = getPoints(edgesToKeep);
    ArrayList<Edge> newEdges = calculKruskal(newPoints);
    for(Edge e : newEdges) {
    	if(e.distance()>edgeThreshold) {
    		System.err.println("err");
    	}
    }
    return newEdges;
  }
  //-----------------------------------------------------------------------------------------------------------------
  //-------------------------------------------- STEINER WITHOUT BUDGET ---------------------------------------------
  //-----------------------------------------------------------------------------------------------------------------

  public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
    int[][] paths = calculShortestPaths(points, edgeThreshold);

    ArrayList<Edge> steiner = steiner(points, edgeThreshold, hitPoints,paths);

    System.out.println("Total distance : " + totalDistance(steiner)); // 2707             //----------------------debug print

    Tree2D steinerTree = edgesToTree(steiner, steiner.get(0).p);
    return steinerTree;
  }

  //-----------------------------------------------------------------------------------------------------------------
  //---------------------------------------------- STEINER WITH BUDGET ----------------------------------------------
  //-----------------------------------------------------------------------------------------------------------------

  public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) throws Exception {
    int[][] paths = calculShortestPaths(points, edgeThreshold);
    ArrayList<Edge> steiner = steiner(points, edgeThreshold, hitPoints, paths);

    this.firstPoint = hitPoints.get(0); // set the house point

    //get the list of points from the list of edges after steiner
    ArrayList<Point> pt = getPoints(steiner);
    // System.out.println("pt " + pt.size());                                    //----------------------debug print     
    
    System.out.println("begin");                                      //----------------------debug print
    // begin optimiser to get rid of all edges that are the most expensive (the most far from the house = first point in hitPoints)
    ArrayList<Edge> newOptimizeEdge = optimiser(points, paths, pt, hitPoints, edgeThreshold, BUDGET);
  
    hitPoints.sort(comp);
    // System.out.println(hitPoints);                                      //----------------------debug print
    // System.out.println(firstPoint); 
    
    System.out.println("Total distance before opti: " + totalDistance(steiner));
    System.out.println("Total distance after opti: " + totalDistance(newOptimizeEdge));

    //check if the result is valid (within the edgeThreshold)
    for(Edge e : steiner) {
    	if(e.distance()>edgeThreshold) {
        throw new Exception ("distance default: " + edgeThreshold + " distance: " + e.distance());
    	}
    }
    
    Tree2D T = edgesToTree(newOptimizeEdge, firstPoint);
    return T;
  }


// HELPER FUNCTIONS -------------------------------------------------------------------------------------------------

  /**
   * Get the matrix of the shortest paths from each point to each other
   * @param points the list of points
   * @param edgeThreshold the maximum distance between two points to create an edge
   * @return int[][] the matrix of the shortest paths
   */
  public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
    int[][] paths=new int[points.size()][points.size()];
    double[][] dist=new double[points.size()][points.size()];

    //initialise paths and dist matrix
    for (int i=0;i<paths.length;i++) for (int j=0;j<paths.length;j++){
      double distance = points.get(i).distance(points.get(j));
      if (distance<edgeThreshold){
        dist[i][j]= distance;
        paths[i][j] = j;
      }else
        dist[i][j]= Double.MAX_VALUE;
    }
    //floyd-warshall
    for (int k=0;k<points.size();k++){
      for (int i=0;i<points.size();i++){
        for (int j=0;j<points.size();j++){
          if (dist[i][j]> dist[i][k] + dist[k][j]){
            dist[i][j]=dist[i][k] + dist[k][j];
            paths[i][j]=paths[i][k] ;
          }
        }
      }
    }
    return paths;
  }

  /**
   * Get the list of points forming the shortest path between two points 
   * @param points the list of points
   * @param paths the matrix of the shortest paths
   * @param i the index of the first point
   * @param j the index of the second point
   * @return ArrayList<Point> the list of points forming the shortest path between i and j
   */ 
  private ArrayList<Point> getPathPoints(ArrayList<Point> points,int[][] paths, int i, int j){
		int k = paths[i][j];
		ArrayList<Point> result = new ArrayList<Point>();
		result.add(points.get(i));
		while (k!=j){
			result.add(points.get(k));
			k = paths[k][j];
		}
		result.add(points.get(j));
		return result;
	}

  /**
   * Get the list of points from the list of edges
   * @param edges the list of edges 
   * @return ArrayList<Point> the list of points
   */
  private ArrayList<Point> getPoints(ArrayList<Edge> edges) {
    ArrayList<Point> points = new ArrayList<Point>();
    for (Edge e : edges) {
      if (!points.contains(e.p))
        points.add(e.p);
      if (!points.contains(e.q))
        points.add(e.q);
    }
    return points;
  }

  /**
   * Optimise the result of the steiner tree
   * This function is used to get rid of the most expensive edges. 
   * Using the backtracking method, it removes the farthest point and tries recursively. 
   * If the removal of this point is not improving, it removes the next one until a limit is reached (23 points to keep).
   * @param points the list of points
   * @param paths the matrix of the shortest paths
   * @param newPoints the list of points after steiner
   * @param hitPoints the list of hitPoints
   * @param edgeThreshold the maximum distance between two points to create an edge
   * @param bidgetLimit the maximum distance of the tree
   * @return ArrayList<Edge> the list of edges of the tree after optimization
   * @throws Exception
   */
  private ArrayList<Edge> optimiser(ArrayList<Point> points, int[][] paths, ArrayList<Point> newPoints, ArrayList<Point> hitPoints, int edgeThreshold, int bidgetLimit) throws Exception {
      ArrayList<Edge> res = calculKruskal(newPoints);
      double total = totalDistance(res);
      int val = checkingS(newPoints, hitPoints);

      //first exit condition : if the number of points is less than MIN_POINTS (23), we stop the recursion
      if (val <= MIN_POINTS) {
          return calculKruskal(newPoints);
      }
      //second exit condition : if the total distance is less than the budget, we stop the recursion
      if (total <= bidgetLimit) {
          return res;
      }
      //test optimisation : print the total distance of the tree after each optimization call
      if (total < maxOpti) {
          maxOpti = total;
          System.out.println("maxOpti                => " + maxOpti);
      }
      
      //sort the list of points by distance to the first point in hitPoints (farthest first)
      hitPoints.sort(comp);

      int range = hitPoints.size();
      ArrayList<Edge> next = null;
      for (int i = 0; i < range && i < 5; i++) {
        //take out the first element of the list which is the farthest point since the sort is in descending order, 
        //then we do the basic processing, calling the steiner function and the optimiser function recursively
          Point p = hitPoints.get(0);
          newPoints.remove(p);
          hitPoints.remove(p);

          ArrayList<Edge> steiner = steiner(points, edgeThreshold, hitPoints, paths);
          
          ArrayList<Point> nextPoints = getPoints(steiner);
          
          next = optimiser(points, paths, nextPoints, hitPoints, edgeThreshold, bidgetLimit);
          double totalTMP = totalDistance(next);
          //if we find a solution we stop the loop 
          if (totalDistance(next) <= bidgetLimit) {
            maxOpti = totalTMP;
            System.out.println("maxOpti-Back                => " + maxOpti);
            // return next;
            break;
          }
          //else we put the point back at the end of the list and we start again
          newPoints.add(p);
          hitPoints.add(p);
      }

      return next;
  }

  /**
   * Count the number of hitPoints in the points
   * This function is used to check if the optimization is valid
   * @param points the list of points
   * @param hitPoints the list of hitPoints
   * @return int the number of hitPoints in the points
   */
  private int checkingS (ArrayList<Point> points, ArrayList<Point> hitPoints) {
    int compt = 0;
    for (Point p : points) {
      if (hitPoints.contains(p)) {
        compt++;
      }
    }
    return compt;
  }

  /**
   * Get the total distance of the list of edges of the tree
   * @param points the list of edges
   * @return double the total distance of the list of edges
   */
  private double totalDistance(ArrayList<Edge> points) {
    double distance = 0;
    for (Edge e : points) {
      distance += e.distance();
    }
    return distance;
  }

  //-----------------------------------------------------------------------------------------------------------------
  //------------------------------------------ KRUSKAL FROM TME 4 ---------------------------------------------------
  //-----------------------------------------------------------------------------------------------------------------

  public ArrayList<Edge> calculKruskal(ArrayList<Point> points) {
    //KRUSKAL ALGORITHM, NOT OPTIMAL FOR STEINER!
    ArrayList<Edge> edges = new ArrayList<Edge>();
    for (Point p: points) {
      for (Point q: points) {
        if (p.equals(q) || contains(edges,p,q)) continue;
        edges.add(new Edge(p,q));
      }
    }
    edges = sort(edges);

    ArrayList<Edge> kruskal = new ArrayList<Edge>();
    Edge current;
    NameTag forest = new NameTag(points);
    while (edges.size()!=0) {
      current = edges.remove(0);
      if (forest.tag(current.p)!=forest.tag(current.q)) {
        kruskal.add(current);
        forest.reTag(forest.tag(current.p),forest.tag(current.q));
      }
    }

    return kruskal;
  }
  private boolean contains(ArrayList<Edge> edges,Point p,Point q){
    for (Edge e:edges){
      if (e.p.equals(p) && e.q.equals(q) ||
          e.p.equals(q) && e.q.equals(p) ) return true;
    }
    return false;
  }
  private Tree2D edgesToTree(ArrayList<Edge> edges, Point root) {
    ArrayList<Edge> remainder = new ArrayList<Edge>();
    ArrayList<Point> subTreeRoots = new ArrayList<Point>();
    Edge current;
    while (edges.size()!=0) {
      current = edges.remove(0);
      if (current.p.equals(root)) {
        subTreeRoots.add(current.q);
      } else {
        if (current.q.equals(root)) {
          subTreeRoots.add(current.p);
        } else {
          remainder.add(current);
        }
      }
    }

    ArrayList<Tree2D> subTrees = new ArrayList<Tree2D>();
    for (Point subTreeRoot: subTreeRoots) subTrees.add(edgesToTree((ArrayList<Edge>)remainder.clone(),subTreeRoot));

    return new Tree2D(root, subTrees);
  }

  private ArrayList<Edge> sort(ArrayList<Edge> edges) {
    if (edges.size()==1) return edges;

    ArrayList<Edge> left = new ArrayList<Edge>();
    ArrayList<Edge> right = new ArrayList<Edge>();
    int n=edges.size();
    for (int i=0;i<n/2;i++) { left.add(edges.remove(0)); }
    while (edges.size()!=0) { right.add(edges.remove(0)); }
    left = sort(left);
    right = sort(right);

    ArrayList<Edge> result = new ArrayList<Edge>();
    while (left.size()!=0 || right.size()!=0) {
      if (left.size()==0) { result.add(right.remove(0)); continue; }
      if (right.size()==0) { result.add(left.remove(0)); continue; }
      if (left.get(0).distance() < right.get(0).distance()) result.add(left.remove(0));
      else result.add(right.remove(0));
    }
    return result;
  }
}

class Edge {
  protected Point p,q;
  protected Edge(Point p,Point q){ this.p=p; this.q=q; }
  protected double distance(){ return p.distance(q); }
}

class NameTag {
  private ArrayList<Point> points;
  private int[] tag;
  protected NameTag(ArrayList<Point> points){
    this.points=(ArrayList<Point>)points.clone();
    tag=new int[points.size()];
    for (int i=0;i<points.size();i++) tag[i]=i;
  }
  protected void reTag(int j, int k){
    for (int i=0;i<tag.length;i++) if (tag[i]==j) tag[i]=k;
  }
  protected int tag(Point p){
    for (int i=0;i<points.size();i++) if (p.equals(points.get(i))) return tag[i];
    return 0xBADC0DE;
  }
}