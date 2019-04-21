/*
 * The class is the node object that the binary tree will consist of. Each node 
 * represents an atom with a unique x,y,z-coordinate. Since the only data 
 * structure being used for this program is a binary tree, only a left and right 
 * child, and single parent, will be stored for each node.
 */

import Jama.*;

public class Node {
    
    //index to keep track of which position the node is within the protein tree
    private int index;
    
    //atom that has x,y,z coordinates when added to the protein tree
    private Atom a;
    
    //tracks if the current node has been checked for feasibility yet
    Boolean checked, stop;
    
    //points to the parent and possible children of the node
    private Node parent, leftChild, rightChild;

    //torsion matrix that will help determine the coordinate position for the 
    //ith atom belonging to the ith node. 
    private double[][] cTorsionMatrix;
    
    //empty constructor
    public Node(){}
    
    //constructor
    public Node(int i, Atom aa, Node p, Boolean b, Boolean st) {
        
        this.checked = b;
        this.stop = st;
        this.index = i;
        this.cTorsionMatrix = new double[4][4];
        
        if(this.index > 1)
            this.parent = p;

        this.a = aa;
        this.leftChild = null;
        this.rightChild = null;
    }
    
    //setter for atom
    public void set_atom(Atom aa)
        {this.a = aa;}
    
    //return coordinate object
    public Atom get_atom()
        {return this.a;}
    
    //setter for checked flag
    public void set_checked(Boolean b)
        {this.checked = b;}
    
    //getter for checked flag
    public Boolean get_checked()
        {return this.checked;}    
    
    //setter for checked flag
    public void set_stop(Boolean s)
        {this.stop = s;}
    
    //getter for checked flag
    public Boolean get_stop()
        {return this.stop;}
    
    public int get_index()
        {return this.index;}
    
    public void set_index(int i)
        {this.index = i;}
    
    //adds a child to the current node
    public void set_leftChild(Node n)
        {this.leftChild = n;}
    
    //adds a child to the current node
    public void set_rightChild(Node n) 
        {this.rightChild = n;}
    
    //setter for torsion matrix for node
    public void set_cTorsionMatrix(double[][] ctm) 
        {this.cTorsionMatrix = ctm;}
    
    //getter for torsion matrix
    public double[][] get_cTorsionMatrix()
        {return this.cTorsionMatrix;}
    
    //prints to screen cumulative torsion matrix
    public void print_cTorsionMatrix() {
    
        Matrix m = new Matrix(this.cTorsionMatrix);
        System.out.println("CUMULATIVE TORSION MATRIX");
        for(int i = 0 ; i < 4; i++){
            for(int k = 0; k < 4; k++)
                System.out.print(String.format("%-7.3f", m.get(i, k)));
            
            System.out.println();
        }
        System.out.println();
    }
    
    public void set_parent(Node n)
        {this.parent = n;}
    
    public Node get_parent()
        {return this.parent;}
    
    //getter for left child
    public Node get_leftChild()
        {return this.leftChild;}
    
    //getter for right child
    public Node get_rightChild()
        {return this.rightChild;}
}