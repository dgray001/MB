package mercedes_benz;
// This class defines a particle in the MC code, but not with LJ parameters

public class Particle {
	
	// default values
	String ID = "MB";
	double x_coordinate = 0.0;
	double y_coordinate = 0.0;
	double orientation = 0.0;
	
	Particle() {
		// no arg constructor
	}
	
	// Other constructor
	Particle(String id, double x, double y, double o) {
		this.ID = id;
		this.x_coordinate = x;
		this.y_coordinate = y;
		this.orientation = o;
	}
	
	public void print_particle() {
		System.out.println("[" + this.getID() + ", " + Double.toString(this.getX()) + ", " + Double.toString(this.getY()) + ", " + Double.toString(this.getO()) + "]\n");
	}
	
	public String getID() {
		return ID;
	}
	public void setID(String iD) {
		this.ID = iD;
	}
	
	public double getX() {
		return x_coordinate;
	}
	public void setX(double x_coordinate) {
		this.x_coordinate = x_coordinate;
	}
	
	public double getY() {
		return y_coordinate;
	}
	public void setY(double y_coordinate) {
		this.y_coordinate = y_coordinate;
	}
	
	public double getO() {
		return orientation;
	}
	public void setO(double orientation) {
		this.orientation = orientation;
	}
}
