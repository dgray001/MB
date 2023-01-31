//package Java_Files;

import java.util.Arrays;

public class MD_utils {
	public static void print_particles_MD(Particle_MD[] AP) {
		System.out.println("[ID, Xcoor, Ycoor, Ocoor, x_vel, y_vel, o_vel, x_for, y_for, o_for, mass, dist.]");
		for (int i = 0; i < AP.length; i++) {
			String[] particle = new String[12];
			particle[0] = AP[i].getID();
			particle[1] = Double.toString(Math.round(AP[i].getX()*1000)/1000.0);
			particle[2] = Double.toString(Math.round(AP[i].getY()*1000)/1000.0);
			particle[3] = Double.toString(Math.round(AP[i].getO()*1000)/1000.0);
			particle[4] = Double.toString(Math.round(AP[i].getVel()[0]*1000/1000.0));
			particle[5] = Double.toString(Math.round(AP[i].getVel()[1]*1000/1000.0));
			particle[6] = Double.toString(Math.round(AP[i].getVel()[2]*1000/1000.0));
			particle[7] = Double.toString(Math.round(AP[i].getFor()[0]*1000/1000.0));
			particle[8] = Double.toString(Math.round(AP[i].getFor()[1]*1000/1000.0));
			particle[9] = Double.toString(Math.round(AP[i].getFor()[2]*1000/1000.0));
			particle[10] = Double.toString(Math.round(AP[i].getM()*1000/1000.0));
			particle[11] = Double.toString(Math.round(AP[i].getD()*1000/1000.0));
			System.out.println(Arrays.toString(particle));
		}
	}
	public static Particle_MD[] setup_particles_MD(int n, int c, int s, double massN, double massC, double massS, double[] q, double[] sigma) {
		Particle_MD[] AP = new Particle_MD[(n+s+c)];
		for (int i = 0; i < (n+s+c); i++) {
			Particle_MD particle = new Particle_MD();
			if (i < c) {
				particle.setID("LJ_c");
				particle.setO(0.0);
				particle.setM(massC);
			}
			else if (i < (c+s)) {
				particle.setID("LJ_s");
				particle.setO(0.0);
				particle.setM(massS);
			}
			else {
				particle.setID("MB");
				particle.setO(Math.random()*Math.PI);
				particle.setM(massN);
			}
			boolean tooClose = true;
			int tries = 0;
			double[] bestTry = {0,0,0};
			while (tooClose) {
				tooClose = false;
				particle.setX(Math.random()*q[0]);
				particle.setY(Math.random()*q[1]);
				for (int j = 0; j < i; j++) {
					double[] u_hat_not = MC_utils.nearest_image(particle, AP[j], q);
					double r_ij = Math.sqrt(u_hat_not[0]*u_hat_not[0] + u_hat_not[1]*u_hat_not[1]);
					double contactDistance = (sigma[i] + sigma[j]) / 2.0;
					if (r_ij < contactDistance) {
						tooClose = true;
						if (r_ij > bestTry[2]) {
							bestTry[0] = particle.getX();
							bestTry[1] = particle.getY();
							bestTry[2] = r_ij;
						}
					}
				}
				tries ++;
				if (tries >= 100000*(n+s+c)) {
					particle.setX(bestTry[0]);
					particle.setY(bestTry[1]);
					System.out.println("Particle too close.");
					break;
				}
			}
			AP[i] = particle;
		}
		return AP;
	}
	public static void Boltzmann_distribution(Particle_MD[] AP, int n, double I, double T) {
		double[][] velocities = new double[AP.length][3];
		for (int i = 0; i < velocities.length; i++) {
			for (int j = 0; j < velocities[0].length; j++) {
				velocities[i][j] = Math.random();
			}
		}
		double[] vSums = {0.0, 0.0, 0.0};
		for (int i = 0; i < velocities.length; i++) {
			if (AP[i].ID == "MB") {
				AP[i].velocity[0] = PPF(velocities[i][0])*Math.sqrt(T/AP[i].mass);
				AP[i].velocity[1] = PPF(velocities[i][1])*Math.sqrt(T/AP[i].mass);
				AP[i].velocity[2] = PPF(velocities[i][2])*Math.sqrt(T/I);
			}
			else {
				AP[i].velocity[0] = PPF(velocities[i][0])*Math.sqrt(T/AP[i].mass);
				AP[i].velocity[1] = PPF(velocities[i][1])*Math.sqrt(T/AP[i].mass);
				AP[i].velocity[2] = 0.0;
			}
			vSums[0] += AP[i].velocity[0];
			vSums[1] += AP[i].velocity[1];
			vSums[2] += AP[i].velocity[2];
		}
		vSums[0] /= AP.length;
		vSums[1] /= AP.length;
		if (n > 0) {
			vSums[2] /= n;
		}
		for (int i = 0; i < AP.length; i++) {
			AP[i].velocity[0] -= vSums[0];
			AP[i].velocity[1] -= vSums[1];
			if (AP[i].ID == "MB") {
				AP[i].velocity[2] -= vSums[2];
			}
		}
	}
	public static double PPF(double x) {
		double result = 1.41421*invErf(2*x-1);
		return result;
	}
	public static double invErf(double x) {
		double result = 0.0;
		double[] coefficients = invErfCo(99);
		for (int i = 0; i < 100; i++) {
			result += (coefficients[i]/(2*i + 1))*Math.pow((0.886227*x),(2*i + 1));
		}
		return result;
	}
	public static double[] invErfCo(int x) {
		double[] coefficients = new double[x+1];
		coefficients[0] = 1.0;
		for (int i = 1; i <= x; i++) {
			double sum = 0.0;
			for (int j = 0; j < i; j++) {
				sum += (coefficients[i]*coefficients[i - 1 - j])/((i + 1)*(2*i + 1));
			}
			coefficients[i] = sum;
		}
		return coefficients;
	}
	
	public static Particle_MD[] velocity(Particle_MD[] AP, double dt, double[][] last_acc, double IHB) { // for full steps
		for (int i = 0; i < AP.length; i++) {
			AP[i].velocity[0] += 0.5*(AP[i].force[0] + last_acc[i][0])*dt/AP[i].mass;
			AP[i].velocity[1] += 0.5*(AP[i].force[1] + last_acc[i][1])*dt/AP[i].mass;
			AP[i].velocity[2] += 0.5*(AP[i].force[2] + last_acc[i][2])*dt/IHB;
		}
		return AP;
	}
	public static Particle_MD[] velocity_half(Particle_MD[] AP, double dt, double IHB) { // for half steps
		for (int i = 0; i < AP.length; i++) {
			AP[i].velocity[0] += 0.5*AP[i].force[0]*dt/AP[i].mass;
			AP[i].velocity[1] += 0.5*AP[i].force[1]*dt/AP[i].mass;
			AP[i].velocity[2] += 0.5*AP[i].force[2]*dt/IHB;
		}
		return AP;
	}
	public static Particle_MD[] position(Particle_MD[] AP, double dt, double[] q, double IHB) {
		for (int i = 0; i < AP.length; i++) {
			double xDist = AP[i].velocity[0]*dt + 0.5*AP[i].force[0]*dt*dt/AP[i].mass; // distance moving in x direction
			double yDist = AP[i].velocity[1]*dt + 0.5*AP[i].force[1]*dt*dt/AP[i].mass;
			AP[i].x_coordinate += xDist; // q = q_not + v*t + .5*a*t^2
			AP[i].y_coordinate += yDist;
			AP[i].orientation += AP[i].velocity[2]*dt + 0.5*AP[i].force[2]*dt*dt/IHB;
			wrap(AP[i], q);
			AP[i].d += xDist*xDist + yDist*yDist; // msd
		}
		return AP;
	}
	public static Particle_MD[] position_fast(Particle_MD[] AP, double dt, double[] q, double IHB) {
		for (int i = 0; i < AP.length; i++) {
			double xDist = AP[i].velocity[0]*dt; // distance moving in x direction
			double yDist = AP[i].velocity[1]*dt;
			AP[i].x_coordinate += xDist; // q = q_not + v*t
			AP[i].y_coordinate += yDist;
			AP[i].orientation += AP[i].velocity[2]*dt;
			wrap(AP[i], q);
			AP[i].d += (xDist*xDist + yDist*yDist); // calculate msd
		}
		return AP;
	}
	public static Particle_MD[] force(Particle_MD[] AP, double[] sigma, double[] epsilon, double[] q, double rHB, double sigHB, double epHB, double Gcon1, double Gcon2) {
		double[] blank = {0.0, 0.0, 0.0};
		for (int i = 0; i < AP.length; i++) {
			AP[i].setFor(blank); // set all forces to 0.0
		}
		for (int i = 0; i < AP.length; i++) {
			for (int j = (i+1); j < AP.length; j++) { // every pair
				double[] u_hat_not = MC_utils.nearest_image(AP[i], AP[j], q);
				double r_ij = Math.sqrt(u_hat_not[0]*u_hat_not[0] + u_hat_not[1]*u_hat_not[1]);
				double sig = 0.5*(sigma[i] + sigma[j]);
				double ep = 0.5*(epsilon[i] + epsilon[j]);
				if (r_ij < 4*sig) {
					double[][] LJf = LJ_force(u_hat_not, r_ij, sig, ep);
					AP[i].setFor(MC_utils.vector_add(AP[i].force, LJf[0]));
					AP[j].setFor(MC_utils.vector_add(AP[j].force, LJf[1]));
				}
				if ((AP[i].ID == "MB")&&(AP[j].ID == "MB")&&(r_ij > 0.5*rHB)&&(r_ij < 1.5*rHB)) {
					double[][] HBf = HB_force(AP[i], AP[j], u_hat_not, r_ij, rHB, sigHB, epHB, Gcon1, Gcon2);
					AP[i].setFor(MC_utils.vector_add(AP[i].force, HBf[0]));
					AP[j].setFor(MC_utils.vector_add(AP[j].force, HBf[1]));
				}
			}
		}
		return AP;
	}
	
	public static void wrap(Particle_MD i, double[] q) { // maintain periodic boundaries
		if (i.x_coordinate > q[0]) {
			i.x_coordinate -= q[0];
		}
		else if (i.x_coordinate < 0.0) {
			i.x_coordinate += q[0];
		}
		if (i.y_coordinate > q[1]) {
			i.y_coordinate -= q[1];
		}
		else if (i.y_coordinate < 0.0) {
			i.y_coordinate += q[1];
		}
	}
	public static double[][] LJ_force(double[] u_not, double rij, double sigma, double epsilon) {
		double pow1 = sigma/rij;
		double pow2 = pow1*pow1;
		double pow6 = pow2*pow2*pow2;
		double force = -24*epsilon*pow6*(2*pow6 - 1)/rij; // dLJ/dr = -24*e*(2*(sigma/r)^12 - (sigma/r)^6)/r
		double fX = force*u_not[0]/rij; // multiply by dr/dx
		double fY = force*u_not[1]/rij; // multiply by dr/dy
		double[][] result = {{fX, fY, 0.0}, {-fX, -fY, 0.0}};
		return result;
	}
	public static double[][] HB_force(Particle_MD i, Particle_MD j, double[] u_not, double rij, double rHB, double sigHB, double epHB, double Gcon1, double Gcon2) {
		double[] u = MC_utils.vector_mult(u_not, 1/rij);
		double[][] i_arms = MC_utils.get_arms(i.orientation);
		double[][] j_arms = MC_utils.get_arms(j.orientation);
		double[] i_dots = {MC_utils.DOT(i_arms[0], u), MC_utils.DOT(i_arms[1], u), MC_utils.DOT(i_arms[2], u)};
		double[] j_dots = {MC_utils.DOT(j_arms[0], u), MC_utils.DOT(j_arms[1], u), MC_utils.DOT(j_arms[2], u)};
		int i_index = 0;
		int j_index = 0;
		for (int ind = 1; ind < 3; ind++) {
			i_index = (i_dots[ind] > i_dots[i_index]) ? ind:i_index;
			j_index = (j_dots[ind] < j_dots[j_index]) ? ind:j_index;
		}
		double i_angle = i.orientation + i_index*2.094395; // angle from unit vector i_hat
		double j_angle = j.orientation + j_index*2.094395;
		double drdx = -u_not[0]/rij; // dr/dx
		double drdy = -u_not[1]/rij; // dr/dy
		double dhidx = -(Math.cos(i_angle) + i_dots[i_index]*drdx)/rij; // dh(i)/dx
		double dhjdx = -(Math.cos(j_angle) + j_dots[j_index]*drdx)/rij; // dh(j)/dx
		double dhidy = -(Math.sin(i_angle) + i_dots[i_index]*drdy)/rij;
		double dhjdy = -(Math.sin(j_angle) + j_dots[j_index]*drdy)/rij;
		double dhi_dphi = (-u_not[0]*Math.sin(i_angle) + u_not[1]*Math.cos(i_angle))/rij; // for torques
		double dhj_dphi = (-u_not[0]*Math.sin(j_angle) + u_not[1]*Math.cos(j_angle))/rij;
		double Gr = MC_utils.Gaussian(rij - rHB, Gcon1);
		double Gi = MC_utils.Gaussian(i_dots[i_index] - 1, Gcon1); // three gaussians
		double Gj = MC_utils.Gaussian(j_dots[j_index] + 1, Gcon1);
		double dGr = (rij - rHB)*Gcon2*Gr;
		double dGi = (i_dots[i_index] - 1)*Gcon2*Gi; // gaussian derivatives
		double dGj = (j_dots[j_index] + 1)*Gcon2*Gj;
		double fX = -epHB*(dGr*drdx*Gi*Gj + Gr*dGi*dhidx*Gj + Gr*Gi*dGj*dhjdx); // Product rule
		double fY = -epHB*(dGr*drdy*Gi*Gj + Gr*dGi*dhidy*Gj + Gr*Gi*dGj*dhjdy);
		double tI = -epHB*Gr*dGi*dhi_dphi*Gj; // torques different
		double tJ = -epHB*Gr*Gi*dGj*dhj_dphi;
		double[][] result = {{fX, fY, tI}, {-fX, -fY, tJ}};
		return result;
	}
	
	public static double virial_pressure(Particle_MD[] AP, double rho, double IHB, double[] q, double[] sigma, double[] epsilon, double rHB, double sigHB, double epHB, double Gcon1, double Gcon2) {
		double pressure = 0;
		double KE = KE_total(AP, IHB);
		for (int i = 0; i < AP.length; i++) {
			for (int j = (i+1); j < AP.length; j++) {
				double[] u_hat_not = MC_utils.nearest_image(AP[i], AP[j], q);
				double r_ij = Math.sqrt(u_hat_not[0]*u_hat_not[0] + u_hat_not[1]*u_hat_not[1]);
				double sig = 0.5*(sigma[i] + sigma[j]);
				double ep = 0.5*(epsilon[i] + epsilon[j]);
				double[][] LJf = LJ_force(u_hat_not, r_ij, sig, ep);
				if ((AP[i].ID == "MB")&&(AP[j].ID == "MB")) {
					double[][] HBf = HB_force(AP[i], AP[j], u_hat_not, r_ij, rHB, sigHB, epHB, Gcon1, Gcon2);
					LJf[0][0] += HBf[0][0];
					LJf[0][1] += HBf[0][1];
				}
				pressure += u_hat_not[0]*LJf[0][0] + u_hat_not[1]*LJf[0][1];
			}
		}
		pressure += 2*KE;
		pressure *= rho/(3*AP.length); // pressure = (rho/3*N)*(2*KE + sum(dot(r_hat, F)))
		return pressure;
	}
	public static double KE_total(Particle_MD[] AP, double IHB) { // find kinetic energy
		double KE = 0.0;
		for (int i = 0; i < AP.length; i++) {
			KE += 0.5*AP[i].mass*AP[i].velocity[0]*AP[i].velocity[0];
			KE += 0.5*AP[i].mass*AP[i].velocity[1]*AP[i].velocity[1];
			KE += 0.5*IHB*AP[i].velocity[2]*AP[i].velocity[2];
		}
		return KE;
	}
	public static void berendsen_thermostat(Particle_MD[] AP, double dt_tau, double kBT, double IHB) {
		double T = 0.66667*KE_total(AP, IHB)/AP.length;
		double scale = Math.sqrt(1 + dt_tau*((kBT/T) - 1)); // scale = sqrt(1 + (dt/tau)*((kBT/Tcurr) - 1))
		for (int i = 0; i < AP.length; i++) {
			AP[i].velocity[0] *= scale;
			AP[i].velocity[1] *= scale;
			AP[i].velocity[2] *= scale;
		}
	}
}
