/*
 * The purpose of this object is to model the distance between two atoms.
 */


public class Distance {
    
    //distance between the two atoms
    private float distance;
    
    //the two atoms that are used to define the distance
    private Atom atom1, atom2;
    
    //empty constructor
    public Distance(){}
    
    //constructor
    public Distance(Atom a1, Atom a2)
        {this.atom1=a1; this.atom2=a2; if(!a1.get_name().equalsIgnoreCase(null))this.calc_distance();}
    
    //METHOD TO CALCULATE//
    //The purpose of this method is to compute the distance between two atoms. 
    //It does so by using the Pythagorean Theorem: 
    //Distance = sqrt[(x1-x2)^2+(y1-y2)^2+(z1-z2)^2].
    private void calc_distance() {
        
        float x = (float)Math.pow((this.atom1.get_xCoord() - this.atom2.get_xCoord()),2);
        float y = (float)Math.pow((this.atom1.get_yCoord() - this.atom2.get_yCoord()),2);
        float z = (float)Math.pow((this.atom1.get_zCoord() - this.atom2.get_zCoord()),2);
    
        //Pythagorean Theorem in R^3
        this.distance = (float)Math.sqrt(x+y+z);
    }
    
    //getter for atom1
    public Atom get_atom1()
        {return this.atom1;}
    
    //getter for atom2
    public Atom get_atom2()
        {return this.atom2;}
    
    //getter for distance
    public float get_distance()
        {return this.distance;}
}