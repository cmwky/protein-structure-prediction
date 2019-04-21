/*
The purpose of this class is to represent an atom, which is described by
* an elemental symbol and an x,y,z-coordinate.
*/

import java.lang.Math.*;

public class Atom {
    
    //x,y,z coordinates for each endpoint (atom) of the bond in R^3
    private double x_coord, y_coord, z_coord;
    
    //elemental symbol for atoms
    private String atom_name;
    
    //index for the ith atom
    private int index;
    
    //constructor
    public Atom (String s, int i, double x, double y, double z) {
        
        //sets elemental symbol for each atom
        this.atom_name = s;

        //sets index
        this.index = i;
        
        //sets x,y,z coords
        this.x_coord = x;
        this.y_coord = y;
        this.z_coord = z;
    }
    
    
    //setter for atom index
    public void set_index(int i)
        {this.index = i;}

    //getter for atom name/symbol
    public String get_name()
        {return this.atom_name;}
    
    //getter for index
    public int get_index()
        {return this.index;}
    
    //getter for x coordinate
    public double get_xCoord()
        {return this.x_coord;}
    
    //getter for y coordinate
    public double get_yCoord()
        {return this.y_coord;}
        
    //getter for z coordinate
    public double get_zCoord()
        {return this.z_coord;}
    
    //setter for x coordinate
    public void set_xCoord(double x)
        {this.x_coord = x;}
    
    //setter for y coordinate
    public void set_yCoord(double y)
        {this.y_coord = y;}
    
    //setter for z coordinate
    public void set_zCoord(double z)
        {this.z_coord = z;}
    
    //setter for all three coordinates
    public void set_xyzCoord(double x, double y, double z) 
        {this.x_coord=x; this.y_coord=y; this.z_coord=z;}
    
    //prints summary of distance info
    public void getATOMinfo() {

        //atom1 information
        System.out.println("Atom #" + this.index + ": " + this.atom_name);
        System.out.println("x: " + String.format("%.10f",this.x_coord) );
        System.out.println("y: " + String.format("%.10f",this.y_coord));
        System.out.println("z: " + String.format("%.10f",this.z_coord));
        System.out.println();
    }
}