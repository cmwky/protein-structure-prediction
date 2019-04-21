/*
 * The purpose of this object is to model the situation where three atoms are
 * being considered to then determine the position of the next, newly appended
 * atom. This collection of three atoms will be considered a "clique", and the
 * atoms will be obtained from the binary tree.
 */


public class Clique {
    
    //angle formed between bonds belonging to atoms i/i+1 and i+1/i+2
    private double theta;
    
    //distances of bonds d_i+1 and d_i+2, where d_i+1 is between atoms i and 
    //i+1, and d_i+2 is between atoms i+1 and i+2, where i is the current atom
    private Distance d1, d2, d3;

    //constructor
    public Clique(int i, Distance dd1, Distance dd2) {
    
        this.d1 = dd1;
        this.d2 = dd2;
        this.d3 = new Distance(this.d1.get_atom1(), this.d2.get_atom2());
        
        if(i==1) 
            this.theta = this.calculate_theta(this.d1.get_distance(), this.d2.get_distance(), this.d3.get_distance());
        if(i==2) 
            this.theta = this.calculate_theta(this.d1.get_distance(), this.d3.get_distance(), this.d2.get_distance());
    }
    
    //getter for theta
    public double get_theta()
        {return this.theta;}
    
    //getter for d1
    public Distance get_distance1()
        {return this.d1;}
    
    //getter for d2
    public Distance get_distance2()
        {return this.d2;}
 
     //getter for d2
    public Distance get_distance3()
        {return this.d3;}   
    
    //the purpose of this method is to take a sequence of three atoms and 
    //determine the angle created between the distances belonging to 
    //atom1&atom2, and atom2 & atom3. 
    //d12 -> distance between atom1 to atom2; d2 -> distance between 
    //atom2 and atom3; d13 -> distance between atom1 and atom3
    private double calculate_theta(double d12, double d23, double d13) {
        
        //law of cosines
        double cos_theta = (Math.pow(d12,2) + Math.pow(d23,2) - Math.pow(d13,2))  
                            / (2*d12*d23);

        return Math.acos(cos_theta);
    }

    //prints to screen the value of theta
    public void print_theta() {
        
        System.out.println("Angle in radians: " + this.theta);
        System.out.println("Angle in degrees: " + Math.toDegrees(this.theta));
    }
}