/*
 * The purpose of this class is to handle all the distance data, beginning with 
 * the coordinates of each atom obtained from the PDB file an to the calculated distances between
 * each pair of atoms 
 */


import java.util.*;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.FileReader;
import java.io.File;

public class DistanceData {

    //stores all atoms from the PDB file
    private ArrayList<Atom> atomsFromFile, onlyBackbone;
   
    //stores the pairwise distances between all pairs of atoms
    private ArrayList<Distance> pairwiseDistances, pairwisePBB;
    
    //distances for covalent bonds within protein backbone
    private ArrayList<Distance> covalentPBB;

    //initializes the ALs containing all atoms read from the file and 
    //their calculated pairwise distance data (simulated NMR).
    public DistanceData() {

        this.atomsFromFile = new ArrayList<>();
        this.onlyBackbone = new ArrayList<>();
        this.pairwisePBB = new ArrayList<>();
        this.covalentPBB = new ArrayList<>();
        this.pairwiseDistances = new ArrayList<>();
    }

    //retrieves the information for each atom from the PDB file
    public void read_atoms_from_file() {

        //creates the path string of the PBD txt file
        String filePath = new File("").getAbsolutePath();
        System.out.println(filePath);
        filePath = filePath.concat("/src/2ncs.txt");
        
        //opens the PDB file
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {

            //used to extract each line (atom) from the PDB file
            String CurrentLine;

            //index assigned to each atom thats extracted from PDB file
            int i = 1;
            
            //reads the PDB file until the end is reached
            while ((CurrentLine = br.readLine()) != null) {

                //checks if the current line being read contains information
                //regarding an atom within the protein
                if (CurrentLine.startsWith("ATOM") /*|| CurrentLine.startsWith("HETATM")*/) {

                    //x,y,z coords for atom
                    float x, y, z;

                    //x,y,z coordinates for atom are read from the 
                    //appropriate columns
                    x = Float.parseFloat(CurrentLine.substring(30, 38));
                    y = Float.parseFloat(CurrentLine.substring(38, 46));
                    z = Float.parseFloat(CurrentLine.substring(46, 54));

                    //a new atom is now created from the recently
                    //obtained information
                    Atom a = new Atom(CurrentLine.substring(12, 16).trim(), i, x, y, z);
                    i++;
                    
                    //newly created atom object is added to the AL of all 
                    //atoms.
                    this.atomsFromFile.add(a);
                }
            }
        } catch (IOException e) 
            {System.out.println("Couldn't load file.");}
    }
    
    //strips the atom list down to only protein backbone atoms
    public void only_proteinBackbone() {
        
        //only stores N, CA and C
        for(int i=0; i < this.atomsFromFile.size()-3; i++)
            if(this.atomsFromFile.get(i).get_name().equalsIgnoreCase("N") &&
               this.atomsFromFile.get(i+1).get_name().equalsIgnoreCase("CA") &&
               this.atomsFromFile.get(i+2).get_name().equalsIgnoreCase("C"))
                for(int k = i; k < i+3; k++)
                    this.onlyBackbone.add(this.atomsFromFile.get(k));  
    }
    
    //takes the atoms belonging to protein backbone and generates an AL
    //that represents the covalent bond distances between each sequential atom.
    public void generate_covalentPBB_distances(ArrayList<Atom> a) {

        for(int i = 0; i < a.size()-1; i++) {
            
            Distance d = new Distance(a.get(i), a.get(i+1));
            this.covalentPBB.add(d);
         }  
    }
    
    //takes the atoms from the PDB file and calculates their pairwise distances
    public ArrayList<Distance> calculate_pairwiseDistances(ArrayList<Atom> a){
    
        ArrayList<Distance> pw = new ArrayList<>();
        
        //beginning with the first atom, calculate the distance between the 
        //current atom and all other atoms. as the index of the current atom
        //increases, only atoms with higher index values will be used to 
        //calculate distances.
        for(Atom aa : a) {
        
            for(int i=a.indexOf(aa)+1; i < a.size(); i++) {
                    
                    //deep copy to a new object
                    Atom a1 = new Atom(aa.get_name(), aa.get_index(), aa.get_xCoord(),aa.get_yCoord(),aa.get_zCoord());
                    a1.set_index(aa.get_index());
                    Atom a2 = new Atom(a.get(i).get_name(), a.get(i).get_index(), a.get(i).get_xCoord(),a.get(i).get_yCoord(),a.get(i).get_zCoord());
                    a2.set_index(a.get(i).get_index());
                    
                    Distance d = new Distance(a1, a2);
                    
                    if(d.get_distance() <= 5.000) {
                        pw.add(d);
                    }
            }
        }    
        return pw;
    }
    
    //getter for protein backbone atoms
    public ArrayList<Atom> get_onlyPBB() 
        {return this.onlyBackbone;}
    
    //setter for covalent bond distances for protein backbone
    public void set_covalentPBB(ArrayList<Distance> cv)
        {this.covalentPBB = cv;}
    
    //getter for covalent bond distances for protein backbone
    public ArrayList<Distance> get_covalentPBB()
        {return this.covalentPBB;}
    
    //getter for pairwise data generated from protein backbone only
    public ArrayList<Distance> get_pairwisePBB()
        {return this.pairwisePBB;}
    
    public void set_pairwisePBB(ArrayList<Distance> pwbb)
        {this.pairwisePBB = pwbb;}
    
    public void set_pairwise_distances(ArrayList<Distance> pw)
        {this.pairwiseDistances = pw;}
    
    //getter for all atoms from PDB file
    public ArrayList<Atom> get_all_atoms()
        {return this.atomsFromFile;}
    
    //getter for all pair-wise distances
    public ArrayList<Distance> get_pairwise_distances()
        {return this.pairwiseDistances;}
    
    //prints to screen all atoms from PDB file
    public void print_allAtoms(){
    
        for(Atom at : this.atomsFromFile)
            at.getATOMinfo();
    }
    
    //prints to screen pairwise generated data based off atoms from PDB file
    public void print_pwDistances() {
    
        for(Distance di : this.pairwiseDistances){
            di.get_atom1().getATOMinfo();
            di.get_atom2().getATOMinfo();
            System.out.println("Distance between atoms: " + di.get_distance());
        }
    }
}