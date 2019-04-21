/*
 * The purpose of this class is to model the currently-being-constructed
 * protein. The data structure that's being used is a binary tree, initialized
 * with three atoms that will have 'fixed' coordinates. 
 */


public class ProteinTree {
    
    //first three nodes contain the first three fixed atoms
    private Node node1, node2, node3;
    
    //empty constructor
    public ProteinTree(){}
    
    //constructor
    public ProteinTree(Clique c)
        {this.set_first3(c);}
    
    //sets the first three nodes in the protein tree (pt), and by doing so,
    //intializing it. the first two nodes only have left children.
    private void set_first3(Clique c) {

        this.create_node1(c);
        this.create_node2(c);
        this.create_node3(c);
        
        this.node1.set_leftChild(this.node2);
        this.node2.set_parent(node1);
        this.node2.set_leftChild(this.node3);
        this.node3.set_parent(node2);
    }
    
    private void create_node1(Clique c){
        
        TorsionMatrix tm = new TorsionMatrix();
        this.node1 = new Node(1, null, null, false, false);
        
        this.node1.set_cTorsionMatrix(tm.generate_B1());
        
        //pulls x,y,z coords from cumulative torsion matrix
        double[] xyz = tm.eval_xyz_fromCTM(this.node1.get_cTorsionMatrix());
        
        Atom a = new Atom(c.get_distance1().get_atom1().get_name(), 
                          1,  
                          
        /*x-coord */    xyz[0],
        /*y-coord */    xyz[1],
        /*z-coord */    xyz[2]);

        //now that the atom has coordinate values, it can be 
        //assigned to the node.
        this.node1.set_atom(a);
 
        this.node1.set_leftChild(null);
        this.node1.set_rightChild(null);
    }
    
    private void create_node2(Clique c){
        
        TorsionMatrix tm = new TorsionMatrix();
        this.node2 = new Node(2, null, null, false, false);
        
        this.node2.set_cTorsionMatrix(tm.eval_cumuTorsionMatrix(this.node1.get_cTorsionMatrix(), tm.generate_B2(c)));
        double[] xyz = tm.eval_xyz_fromCTM(this.node2.get_cTorsionMatrix());
        
        Atom a = new Atom(c.get_distance1().get_atom2().get_name(), 
                          2, 
                          
        /*x-coord */    xyz[0], 
        /*y-coord */    xyz[1], 
        /*z-coord */    xyz[2]);
        
        //now that atom has coordinate values it can be assigned to the node
        this.node2.set_atom(a);
        
        this.node2.set_leftChild(null);
        this.node2.set_rightChild(null);
    }
        
    private void create_node3(Clique c){
        
        TorsionMatrix tm = new TorsionMatrix();
        this.node3 = new Node(3, null, null, false, false);  
        
        this.node3.set_cTorsionMatrix(tm.eval_cumuTorsionMatrix(this.node2.get_cTorsionMatrix(), tm.generate_B3(c)));
        double[] xyz = tm.eval_xyz_fromCTM(this.node3.get_cTorsionMatrix());
        
        Atom a = new Atom(c.get_distance2().get_atom2().get_name(), 
                          3, 
                
        /*x-coord */      xyz[0], 
        /*y-coord */      xyz[1], 
        /*z-coord */      xyz[2]);

        //now that atom has coordinate values, it can be assigned to the node
        this.node3.set_atom(a);
        
        this.node3.set_leftChild(null);
        this.node3.set_rightChild(null);
    }
    
    //sets the ith node. this entails actually creating two nodes, one for each 
    //potentially feasible position for the newly added atom, and making them 
    //the children of the current node. each node will be tested and then 
    //pruned if needed.
    public Node create_nodei(Node curr, int index, Clique c1, Clique c2) {
    
        TorsionMatrix tm = new TorsionMatrix();
        Node n = new Node(curr.get_index()+1, null, curr, false, false);
        
        //assigns cumulative torsion matrix by taking the parents ctm and 
        //mutliplying it with the new nodes torsion matrix.
        n.set_cTorsionMatrix(tm.eval_cumuTorsionMatrix(n.get_parent().get_cTorsionMatrix(), tm.generate_Bi(c1, c2)));
        
        //obtains x,y,z coords from ctm
        double[] xyz = tm.eval_xyz_fromCTM(n.get_cTorsionMatrix());
        
        Atom a = new Atom(c2.get_distance2().get_atom2().get_name(), 
                          index,
                
        /*x-coord */      xyz[0], 
        /*y-coord */      xyz[1], 
        /*z-coord */      xyz[2]);

        //assuming the right node is being created, negate the sign of the 
        //z coordinate for the pre existing left node.
        if(curr.get_leftChild() != null)
               a.set_zCoord(-a.get_zCoord());
        
        //now that the atom is created, it can be assigned to the new node
        n.set_atom(a);
   
        return n;
    }
    
    //getter for node1
    public Node get_node1()
        {return this.node1;}
    
    //getter for node2
    public Node get_node2()
        {return this.node2;}
        
    //getter for node3
    public Node get_node3()
        {return this.node3;}
    
    //prints to screen the first three nodes that initialize the protein tree
    public void print_first3() {
    
        this.node1.get_atom().getATOMinfo();
        this.node2.get_atom().getATOMinfo();
        this.node3.get_atom().getATOMinfo();
    }
}