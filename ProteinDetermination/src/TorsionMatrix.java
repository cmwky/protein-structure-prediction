/*
 * The purpose of this class is to generate the appropriate torsion martix
 * for the ith node added to the currently being constructed protein.
 */

import Jama.*;

public class TorsionMatrix {
    
    public double[][] generate_B1() {
        
        double[][] B1 = new double[][]{ {1,0,0,0},     
                                        {0,1,0,0}, 
                                        {0,0,1,0}, 
                                        {0,0,0,1} }; 
        return B1;
    }
        
    //uses the "first" clique to create torsion matrix for second node.
    public double[][] generate_B2(Clique c){
    
        double[][] B2 = new double[][]{ {-1,0,0, -c.get_distance1().get_distance()},
                                        {0,1,0,0},
                                        {0,0,-1,0},
                                        {0,0,0,1} };
        return B2;
    }
    
    //torsion matrix for third node using "first" clique.
    public double[][] generate_B3(Clique c){
        
        double cos_theta, sin_theta;
        
        cos_theta = Math.cos(c.get_theta());
        sin_theta = Math.sin(c.get_theta());
        
        double[][] B3 = new double[][]{ {-cos_theta, -sin_theta, 0, -c.get_distance2().get_distance()*cos_theta},
                                        {sin_theta, -cos_theta, 0, c.get_distance2().get_distance()*sin_theta},
                                        {0,0,1,0},
                                        {0,0,0,1} };
        
        return B3;
    }
    
    //uses the cliques comprised of atoms i-3 to i to create torsion matrix
    //for the ith atom.
    public double[][] generate_Bi(Clique c1, Clique c2) {
    
        double cos_theta, sin_theta, omega, sin_omega, cos_omega;
        
        omega = this.calculate_omega(c1, c2);
        
        cos_theta = Math.cos(c2.get_theta());
        sin_theta = Math.sin(c2.get_theta());
        
        cos_omega = Math.cos(omega);
        sin_omega = Math.sin(omega);
        
        double[][] Bi = new double[][]{ {-cos_theta, -sin_theta, 0, -c2.get_distance2().get_distance()*cos_theta}, 
                                        {sin_theta*cos_omega, -cos_theta*cos_omega, -sin_omega, c2.get_distance2().get_distance()*sin_theta*cos_omega}, 
                                        {sin_theta*sin_omega, -cos_theta*sin_omega, cos_omega, c2.get_distance2().get_distance()*sin_theta*sin_omega}, 
                                        {0,0,0,1} };
        
        Matrix m = new Matrix(Bi);
        
        m.norm2();
        
        return Bi;
    }
    
    //takes two cliques that have overlapping bonds, and calculates the 
    //angle between the two planes that are created from each of the cliques
    public double calculate_omega(Clique c1, Clique c2) {
    
        //non-bond distances
        double d13, d14, d24;
        
        //cosine of the angle between the two planes that are defined by each
        //clique.
        double cos_omega;

        Distance ddd = new Distance(c1.get_distance1().get_atom1(),c2.get_distance2().get_atom2());
        
        d13 = c1.get_distance3().get_distance();
        d14 = ddd.get_distance();
        d24 = c2.get_distance3().get_distance();
        
        //cosine of the angle between the two planes that is created from 
        //the two overlapping cliques.
        cos_omega = (Math.pow(c1.get_distance1().get_distance(),2) + 
                     Math.pow(d24, 2) - 
                     2 * c1.get_distance1().get_distance() * d24 * Math.cos(c1.get_theta()) * Math.cos(c2.get_theta()) -
                     Math.pow(d14,2))
                    / 
                    (2 * c1.get_distance1().get_distance() * d24 * Math.sin(c1.get_theta() * Math.sin(c2.get_theta())));
        
        //block of code to print variable details, used for verifying by hand
//        System.out.println("r_i+1:  "+ c1.get_distance1().get_distance());
//        System.out.println("r_i+2:  "+ c1.get_distance2().get_distance());
//        System.out.println("r_i+3:  "+ c2.get_distance2().get_distance());
//        System.out.println("d_i_i+2:  "+ d13);
//        System.out.println("d_i_i+3:  "+ d14);
//        System.out.println("d_i+1_i+3:  "+ d24);
//        System.out.println("cosine of theta_i+2:  " + Math.cos(c1.get_theta()));
//        System.out.println("cosine of theta_i+3:  " + Math.cos(c2.get_theta()));
//        System.out.println("sine of theta_i+2:  " + Math.sin(c1.get_theta()));
//        System.out.println("sine of theta_i+3:  " + Math.sin(c2.get_theta()));
//        System.out.println("cos omega:  " + cos_omega);
//        System.out.println("omega in degs:  " + Math.toDegrees(Math.acos(cos_omega)));
//        System.out.println("omega in rads:  " + Math.toRadians(Math.acos(cos_omega)));
//        System.out.println();
        
        return Math.acos(cos_omega);
    }
    
    //this method takes the cumulative torsion matrix from the previous node 
    //and multiplies it with the current nodes torsion matrix to obtain the 
    //cumulative torsion matrix for the current node.
    public double[][] eval_cumuTorsionMatrix(double[][] B, double[][] C) {
    
        //prepares matrices for multiplying
        Matrix m = new Matrix(B);
        Matrix n = new Matrix(C);

        return m.times(n).getArray();
    }
    
    //the purpose of this method is to take the cumulative torsion matrix
    //of a given node and compute the x,y,z coordinates for the atom belonging
    //to that node. 
    public double[] eval_xyz_fromCTM(double[][] ctm) {
    
        //used to return the x,y,z coordinates
        double[] dd = new double[3];
        
        //multiplied with the ctm
        double[][] mat = {{0},{0},{0},{1}};
        
        //preparing matrices to be multiplied
        Matrix m = new Matrix(ctm);
        Matrix n = new Matrix(mat);
        
        //result of the matrix multiplication
        Matrix mn = m.times(n);
        
        //assigns the x,y,z coords to the returning array
        dd[0] = mn.get(0, 0);
        dd[1] = mn.get(1, 0);
        dd[2] = mn.get(2, 0);
        
        return dd;
    }
}