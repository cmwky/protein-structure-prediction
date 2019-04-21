/*
 * The purpose of this class is to take a tree and prune the nodes that contain
 * atoms whose calculated coordinates don't produce distances that agree with
 * the distances obtained from NMR. 
 */

import java.util.*;
//import Jama.*;

public class PruneTree {
   
    //AL that will contain all pairwise distances from the current solution
    private ArrayList<Distance> solutionDistances, NMRDistances;
    private Node currNode;
    double DME;
    
    //empty constructor
    public PruneTree(){}
    
    //constructor
    public PruneTree(Node n, ArrayList<Distance> nmr){
        
        DME = 0;
        this.currNode = n;
        this.NMRDistances = nmr;
        this.solutionDistances = new ArrayList<>();
    }
    
    public void calculate_DME() {
    
        for(Distance sd : this.solutionDistances) {
        
            for(Distance nd : this.NMRDistances) {
            
                if((sd.get_atom1().get_index()==nd.get_atom2().get_index() && sd.get_atom2().get_index()==nd.get_atom1().get_index())) {
                    this.DME += Math.pow(sd.get_distance() - nd.get_distance(), 2);
                }
            }
        }
        
        this.DME = Math.sqrt(this.DME)/(currNode.get_index()-1);
    }
    
    //will take the current node and traverse back to the first atom, with
    //each iteration taking the parent node and placing it into an array, to 
    //then generate pairwise distance data.
    public void generate_solutionPWDistances() {
    
        //same size as parent level for current node
        Node[] solNodes = new Node[this.currNode.get_index()-1];
        
        int i = 0;
        Node n = currNode.get_parent();
        
        //gets all nodes of current solution and adds to solNodes[]
        while(n != null){
            
            solNodes[i] = n;
            n = n.get_parent();
            i++;
        }
        
        //calculates and stores distance between current node's atom and every 
        //other node's atom within the current solution
        for(i = 0; i < solNodes.length; i++) {
            
            Distance d;
            
            if(this.currNode.get_index() - i+1 > 3) {
                
                d = new Distance(this.currNode.get_atom(), solNodes[i].get_atom());
                this.solutionDistances.add(d);
            }
        }
    }
    
    //getter for DME
    public double get_DME()
        {return this.DME;}
}