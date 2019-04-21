/* Charles M Whitehead
 * 
 * The purpose of this program is to take sparse NMR distance data and calculate
 * the 3 dimensional orientation (structure) of the protein that the distance
 * data represents. This is acheived using a "branch and prune" AI technique
 * will adopting specific chemical and physical assumptions.
 */

import java.util.*;
import java.io.*;

public class ProteinDetermination {
    
    //contains all solutions
    ArrayList<Node[]> solutions = new ArrayList<>();
    
    public static void main(String[] args) {
        
        ProteinDetermination pd = new ProteinDetermination();
        
        //this block of code simulates NMR results by taking a PDB file and 
        //extracting the coordinates for each atom, and then calculating the 
        //distance from every atom to every other atom, given the constraint that
        //each distance can not exceed 5 Angstrom, essentially creating 
        //a set of sparse pairwise data. the results of this block of code is an
        //AL of distance objects that contain a distance between two atoms.
        DistanceData d = new DistanceData();
        
        //reads PDB file and creates an AL of atoms
        d.read_atoms_from_file();
        
        //takes the AL of atoms from PDB file and copies only protein backbone
        d.only_proteinBackbone();
        
        //computes pairwise distance data for all atoms
        d.set_pairwise_distances(d.calculate_pairwiseDistances(d.get_all_atoms()));
        
        //computes pairwise data for protein backbone only
        d.set_pairwisePBB(d.calculate_pairwiseDistances(d.get_onlyPBB()));
        
        //takes protein backbone atoms and creates an AL of distances to 
        //represent covalent bonds between sequential atoms
        d.generate_covalentPBB_distances(d.get_onlyPBB());
        
        //constructs the protein backbone
        pd.construct_proteinBB(d.get_covalentPBB(), d.get_covalentPBB());
    }
    
    //the purpose of this method is to construct the backbone of the protein. 
    //this method will later be used to contruct the backbone for each amino 
    //acid side chain.
    private void construct_proteinBB(ArrayList<Distance> covalentDistances, ArrayList<Distance> backbonePWDistances) {
        
        ProteinDetermination pd = new ProteinDetermination();
        
        //used to initialize the backbone, using the first two distances 
        //that is comprised of the first three atoms.
        Clique c = new Clique(1, covalentDistances.get(0), covalentDistances.get(1));
        
        //initializes the protein tree now that the first clique is created,
        //which includes fixing the first three atoms spatial orientation.
        ProteinTree pt = new ProteinTree(c);
        pt.print_first3();
        
        //for each "moment", assign an atom to a new node and then append that
        //node on the backbone tree. the i index represents the number of atoms
        //that are added to the backbone, the k index represents the kth 
        //potential position for that atom.
        for(int i = 0; i < /*covalentDistances.size()*/17; i++) {
            for(int k = 0; k < Math.pow(2, i+1); k++)
                pd.place_atom(pt.get_node3(), covalentDistances, backbonePWDistances, i+3);  
          
            System.out.println("Atom #"+ (i+4) + " is now being placed...");
        }
        
       System.out.println();
       generate_solutions(pt.get_node1(), backbonePWDistances, 20/*covalentDistances.size()*/);
       this.print_solutions();
    }
    
    //prints all solutions from solutions AL
    public void print_solutions() {

        int k = 0;
        for(Node[] sol : this.solutions) {
            
                System.out.println(k);

                try {
                    try (PrintWriter writer = new PrintWriter("solution"+Integer.toString(k++)+".pdb", "UTF-8")) {
                        
                        for(int i = sol.length-1; i >= 0; i--) {

                            //my failed attempt at formatting similiar to that
                            //of a PDB file
                            writer.println("ATOM" + 
                                            String.format("%1$6s", sol[i].get_atom().get_index())+ 
                                                String.format("%1$4s",sol[i].get_atom().get_name()) + "                 " +
                                                String.format("%.4f  ", sol[i].get_atom().get_xCoord()) + 
                                                    String.format("%.4f  ", sol[i].get_atom().get_yCoord()) +
                                                        String.format("%.4f", sol[i].get_atom().get_zCoord()));
                        }
                        writer.close();
                    }
                } 
                catch (IOException e) {}
        }
    }
    
    //generates all solutions
    public void generate_solutions(Node n, ArrayList<Distance> nmr, int i) {
    
        PruneTree pt;
        
        if(n.get_index()==i) {
        
            pt = new PruneTree(n, nmr);
            pt.calculate_DME();
            System.out.println((double)pt.get_DME());
            Node nn = n;
            Node[] sol = new Node[i];
            
            int k = 0;
            while(nn != null) { 
                
                sol[k] = nn;
                nn = nn.get_parent();
                k++;
            }
            this.solutions.add(sol);
        }
        else { 
        
            if(n.get_leftChild() != null && !n.get_leftChild().get_stop())
                this.generate_solutions(n.get_leftChild(), nmr, i);
            
            if(n.get_rightChild() != null && !n.get_rightChild().get_stop())
                this.generate_solutions(n.get_rightChild(), nmr, i);
            }     
        }
    
    //searches for a node that is the same level as the should-be parent, one
    //that also doesn't already have children, and adds the new node (atom)
    //as a child to that node.
    public void place_atom(Node currNode, ArrayList<Distance> covDistances, ArrayList<Distance> backbonePWDistances, int i) {
        
        ProteinTree pt = new ProteinTree();

        if(currNode.get_index() < i && !currNode.get_stop()) {

            if(!currNode.get_leftChild().get_stop())
                place_atom(currNode.get_leftChild(), covDistances, backbonePWDistances, i);
            
            if(currNode.get_rightChild() != null && !currNode.get_rightChild().get_stop())
                place_atom(currNode.get_rightChild(), covDistances, backbonePWDistances, i);
            
        }
        else
        if(currNode.get_index() == i && !currNode.get_checked() && !currNode.get_stop()) {
            
            if(currNode.get_rightChild()==null) {

                Clique[] moment = this.create_moment(currNode, covDistances.get(i-1));
                
                Node newNode = pt.create_nodei(currNode, covDistances.get(i).get_atom1().get_index(), moment[0], moment[1]);

                if(currNode.get_leftChild()==null) {

                    currNode.set_leftChild(newNode);
                    System.out.println("::NEWLY ADDED ATOM::");
                    System.out.println("Added to the left!");
                    System.out.println("parentNodeLevel: "+ currNode.get_index());
                    System.out.println("parentNodeAtom: " + currNode.get_atom().get_name());
                    newNode.get_atom().getATOMinfo();
                    newNode.print_cTorsionMatrix();
                }
                else {
                    
                    currNode.set_rightChild(newNode); 
                   System.out.println("::NEWLY ADDED ATOM::");
                    System.out.println("Added to the right!");
                    System.out.println("parentNodeLevel: "+ currNode.get_index());
                   System.out.println("parentNodeAtom: " + currNode.get_atom().get_name());
                   newNode.get_atom().getATOMinfo();
                    newNode.print_cTorsionMatrix();
                }                   
            }
            
            //if the current node has no children
            if(currNode.get_leftChild() != null && 
                        currNode.get_rightChild() != null) {
                
                double dme1, dme2;
                PruneTree prune1, prune2;
                double epsilon = .001;
                
                //left child DME
                prune1 = new PruneTree(currNode.get_leftChild(), backbonePWDistances);
                prune1.generate_solutionPWDistances();
                prune1.calculate_DME();
                dme1 = prune1.get_DME();
                
                //right child DME
                prune2 = new PruneTree(currNode.get_rightChild(), backbonePWDistances);
                prune2.generate_solutionPWDistances();
                prune2.calculate_DME();
                dme2 = prune2.get_DME();
                
                //prunes the left branch
                if(dme1 > epsilon) {
                    currNode.get_leftChild().set_stop(true);
                    System.out.println("DME1: "+ dme1);
                    System.out.println("left branch trimmed");
                }
                
                //prunes the right branch
                if (dme2 > epsilon) {
                    currNode.get_rightChild().set_stop(true);
                    System.out.println("DME2: "+ dme2);
                    System.out.println("right branch trimmed");
                }
                
                currNode.set_checked(true);
            }
        }
    }
    
    //the purpose of this method is to take a distance and 
    public Distance displace_distance(Distance d, Atom a) {
        
        //stores the difference in each component
        double xDiff, yDiff, zDiff;
        
        //gets the differences in component values so that atom2 can be shifted
        //equally as atom1.
        xDiff = d.get_atom1().get_xCoord() - a.get_xCoord();
        yDiff = d.get_atom1().get_yCoord() - a.get_yCoord();
        zDiff = d.get_atom1().get_zCoord() - a.get_zCoord();
        
        //distance is attached to a, so to speak
        d.get_atom1().set_xyzCoord(a.get_xCoord(),a.get_yCoord(), a.get_zCoord());
        
        //atom is updated accordingly so that the value of the 
        //distance is maintained.
        d.get_atom2().set_xyzCoord(d.get_atom2().get_xCoord() - xDiff,
                                    d.get_atom2().get_yCoord() - yDiff, 
                                    d.get_atom2().get_zCoord() - zDiff);
        return d;
    }
    
    //this method takes the position in which the newly added atom is to be, 
    //and gets the positions of the last three parents of the current node. 
    public Clique[] create_moment(Node currNode, Distance currDist){
            
        //will store the last three parents of the node currently being added
        Node[] parentNodes = new Node[3];
        
        //two cliques that will be used to construct the moment
        Clique c1,c2;
        
        //d1 and d2 come from the tree structure, d3 comes from the covalent
        //bond distances AL
        Distance d1, d2, d3;
        
        //array to return the two cliques of the moment
        Clique[] moments = new Clique[2];
        
        //used for deep copy
        Atom a1, a2;
        
        //deep copy
        a1 = new Atom(currDist.get_atom1().get_name(),
                        currDist.get_atom1().get_index(),
                        currDist.get_atom1().get_xCoord(),
                        currDist.get_atom1().get_yCoord(),
                        currDist.get_atom1().get_zCoord());
        
        //deep copy
        a2 = new Atom(currDist.get_atom2().get_name(), 
                        currDist.get_atom2().get_index(),
                        currDist.get_atom2().get_xCoord(),
                        currDist.get_atom2().get_yCoord(),
                        currDist.get_atom2().get_zCoord());
        
        //deep copy of distance defined by a1 and a2
        Distance updatedDist = new Distance(a1, a2);

        //assigns the last three parents to array, easier to handle
        parentNodes[0] = currNode.get_parent().get_parent();  
        parentNodes[1] = currNode.get_parent();       
        parentNodes[2] = currNode;
       
        //creates new distances from the tree structure atoms
        d1 = new Distance(parentNodes[0].get_atom(), parentNodes[1].get_atom());
        d2 = new Distance(parentNodes[1].get_atom(), parentNodes[2].get_atom());
        
        //displaces the distance between the parent of the newly added atom 
        //and the current atom.
        this.displace_distance(updatedDist, parentNodes[2].get_atom());
        
        //create two cliques using atoms i-3 -> i, with the ith atom being 
        //added.
        c1 = new Clique(1,d1,d2);
        c2 = new Clique(2,d2,updatedDist);
        
        //add cliques to moment array to then return
        moments[0] = c1;
        moments[1] = c2;
        
        return moments;
    }
}