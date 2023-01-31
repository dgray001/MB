//package Java_Files;
// This class extends the MC particle class to add velocity, force, mass, and distance traveled

public class Particle_MD extends Particle {

	// Default values
	double[] velocity = {0,0,0};
	double[] force = {0,0,0};
	double mass = 1.0;
	double d = 0.0; // distance traveled
	
	// Constructors
	Particle_MD() {
	}
	Particle_MD(double[] vel, double[] acc, double m) {
		this.velocity = vel;
		this.force = acc;
		this.mass = m;
	}
	
	// Getters and Setters
	public double[] getVel() {
		return velocity;
	}
	public void setVel(double[] vel) {
		this.velocity = vel;
	}
	public double[] getFor() {
		return force;
	}
	public void setFor(double[] acc) {
		this.force = acc;
	}
	public double getM() {
		return mass;
	}
	public void setM(double m) {
		this.mass = m;
	}
	public double getD() {
		return d;
	}
	public void setD(double d) {
		this.d = d;
	}
}
