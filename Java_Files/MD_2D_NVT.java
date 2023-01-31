package mercedes_benz;

/* Daniel Gray
 * Goal: recreate python code in java of MD MB water to increase speed
 * 
 * Change Log:
 * 12-July-2018:
 * Created File
 * 16-July-2018
 * Made outline of program
 * Made some functions in MD_utils
 * 17-July-2018
 * Finished first draft of program
 * 
 * Pseudocode:
 * Parameters
 * Initialization
 * Main Loop
 * Data Analysis and Output
 */

import java.util.ArrayList;
import java.util.Arrays;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

public class MD_2D_NVT {

	public static void main(String[] args) {
		
		int n = 0; // MB waters
		int c = 0; // large LJ
		int s = 0; // small LJ
		double time_total = 0; // time
		double rho = 1.0; // density
		double kBT = 0.0; // temperature
		boolean NVE = false; // turns of thermostat
//		boolean restart = false; // not functional right now
		n = MC_utils.choose(n,  "Enter # of MB waters:");
		c = MC_utils.choose(c, "Enter # of large LJ particles:");
		s = MC_utils.choose(s, "Enter # of small LJ particles:");
		kBT = MC_utils.choose(kBT, "Enter the temperature:");
		time_total = MC_utils.choose(time_total, "Enter the simulation time:");
		rho = MC_utils.choose(rho, "Enter the density:");
		NVE = MC_utils.choose(NVE, "Press 'y' for NVE (anything else to not)");
//		restart = MC_utils.choose(restart, "Enter 'y' for restart (anything else to not)");
		
		final long start_time = System.currentTimeMillis();	
		
		// Parameters
		final int N = n + c + s;
		final double time_step = 0.001; // dt
		double time_current = 0.0; // current time
		int currentStep = 1; // to have an int
		final double[] q = {Math.sqrt(N/rho), Math.sqrt(N/rho)}; // box
		
		final String filePath = "C:/Users/Yo-Dan/Documents/School/Internships/Summer 2018-Duquesne/Java_MD_Data/";
		
		final double mass_n = 1.0; // 18 g/mol
		final double sigma_n = 0.7;
		final double epsilon_n = 0.1; // MB water
		final double mass_c = 4900; // 88,200 g/mol
		final double sigma_c = 49.0;
		final double epsilon_c = 0.1; // LJ_C parameters
		final double mass_s = 1.0;
		final double sigma_s = 0.7;
		final double epsilon_s = 0.1; // LJ_S parameters
		
		final double r_HB = 1.0; // ideal H_bond distance
		final double epsilon_HB = -1.0; // H_bond strength
		final double sigma_HB = 0.085; // H_bond spread
		final double G_constant1 = -1/(2*sigma_HB*sigma_HB); // to speed up calculations
		final double G_constant2 = -1/(sigma_HB*sigma_HB);
		final double I_HB = 0.0126; // MB particle moment of inertia
		final double tau = 0.20; // coupling constant for Berendsen thermostat
		
		ArrayList<Double> PE_list = new ArrayList<> (); // data lists
		ArrayList<Double> KE_list = new ArrayList<> ();
		ArrayList<Double> U_list = new ArrayList<> ();
		ArrayList<Double> temp_list = new ArrayList<> ();
		ArrayList<Double> pres_list = new ArrayList<> ();
		ArrayList<Double> heat_list = new ArrayList<> ();
		ArrayList<double[]> rdf_list = new ArrayList<> ();
		ArrayList<Double> msdn_list = new ArrayList<> ();
		ArrayList<Double> msdc_list = new ArrayList<> ();
		ArrayList<Double> msds_list = new ArrayList<> ();
		
		int thermostatStep = 5;
		if (NVE) {
			thermostatStep = (int)(time_total/time_step); // no thermostat in NVE
		}
		final int checkStep = (int)(0.5/time_step); // every 0.5 time units
		int eqStep1 = (int)(100/time_step); // energy and temperature
		int eqStep2 = (int)(900/time_step); // the rest
		if (time_total < eqStep1) {
			eqStep1 = (int)(0.5*time_total/time_step); // automatically adjust equilibration
		}
		else if (time_total < eqStep2) {
			eqStep2 = eqStep1;
		}
		final int printStep = (int)(10/time_step);
		
		double[] sigma = new double[N]; // make LJ parameter arrays
		double[] epsilon = new double[N];
		for (int i = 0; i < N; i++) {
			if (i < c) {
				sigma[i] = sigma_c;
				epsilon[i] = epsilon_c;
			}
			else if (i < (c+s)) {
				sigma[i] = sigma_s;
				epsilon[i] = epsilon_s;
			}
			else {
				sigma[i] = sigma_n;
				epsilon[i] = epsilon_n;
			}
		}
		
		// Initialization
		System.out.println("Initiating Particles");
		Particle_MD[] ALL_PARTICLES = MD_utils.setup_particles_MD(n, c, s, mass_n, mass_c, mass_s, q, sigma);
		MD_utils.Boltzmann_distribution(ALL_PARTICLES, n, I_HB, kBT);
		
		MD_utils.force(ALL_PARTICLES, sigma, epsilon, q, r_HB, sigma_HB, epsilon_HB, G_constant1, G_constant2);
		double PE_current = MC_utils.U_total(ALL_PARTICLES, q, sigma, epsilon, G_constant1, r_HB, epsilon_HB);
		double KE_current = MD_utils.KE_total(ALL_PARTICLES, I_HB);
		MD_utils.print_particles_MD(ALL_PARTICLES);
		System.out.println("Equilibrating. Ui = " + PE_current + " (Ti = " + KE_current*0.6667/N + ")");
		
		PrintWriter data = null;
		try {
			data = new PrintWriter(filePath + "MD_data" + Double.toString(kBT) + ".txt", "UTF-8");
		}
		catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		
		// Main Loop
		while (time_current <= time_total) {
			ALL_PARTICLES = MD_utils.velocity_half(ALL_PARTICLES, time_step, I_HB);
			ALL_PARTICLES = MD_utils.position(ALL_PARTICLES, time_step, q, I_HB);
			ALL_PARTICLES = MD_utils.velocity_half(ALL_PARTICLES, time_step, I_HB);
			ALL_PARTICLES = MD_utils.force(ALL_PARTICLES, sigma, epsilon, q, r_HB, sigma_HB, epsilon_HB, G_constant1, G_constant2);
			
			if (currentStep%thermostatStep == 0) {
				MD_utils.berendsen_thermostat(ALL_PARTICLES, time_step/tau, kBT, I_HB);
			}
			
			if (currentStep%checkStep == 0) {
				if (currentStep%printStep == 0) {
					System.out.println("Time: " + Double.toString(Math.round(time_current)));
				}
				if (currentStep > eqStep2) {
					PE_current = MC_utils.U_total(ALL_PARTICLES, q, sigma, epsilon, G_constant1, r_HB, epsilon_HB);
					KE_current = MD_utils.KE_total(ALL_PARTICLES, I_HB);
					PE_list.add(PE_current/N);
					KE_list.add(KE_current/N);
					U_list.add((PE_current + KE_current)/N);
					temp_list.add(0.66667*KE_current/N);
					double pres_current = MD_utils.virial_pressure(ALL_PARTICLES, rho, I_HB, q, sigma, epsilon, r_HB, sigma_HB, epsilon_HB, G_constant1, G_constant2);
					double heat_current = MC_utils.C_p(PE_list, N, kBT);
					double msdn = 0;
					double msdc = 0;
					double msds = 0;
					for (int i = 0; i < N; i++) {
						if (ALL_PARTICLES[i].ID.equals("MB")) {							
							msdn += ALL_PARTICLES[i].d;
						}
						else if (ALL_PARTICLES[i].ID.equals("LJ_c")) {
							msdc += ALL_PARTICLES[i].d;
						}
						else if (ALL_PARTICLES[i].ID.equals("LJ_s")) {
							msds += ALL_PARTICLES[i].d;
						}
					}
					pres_list.add(pres_current);
					heat_list.add(heat_current);
					rdf_list.add(MC_utils.radDist(ALL_PARTICLES, 5*N, Math.min(q[0]/2.0, q[1]/2.0), q));
					if (n > 0) {
						msdn_list.add(msdn/n);
					}
					if (c > 0) {
						msdc_list.add(msdc/c);
					}
					if (s > 0) {
						msds_list.add(msds/s);
					}
					data.println(time_current + " " + PE_current + " " + KE_current + " " + (PE_current+KE_current) + " " + (0.66667*KE_current/N) + " " + pres_current + " " + heat_current + " " + msdn/n + " " + msdc/c + " " + msds);
				}
				else if (currentStep > eqStep1) {
					PE_current = MC_utils.U_total(ALL_PARTICLES, q, sigma, epsilon, G_constant1, r_HB, epsilon_HB);
					KE_current = MD_utils.KE_total(ALL_PARTICLES, I_HB);
					PE_list.add(PE_current/N);
					KE_list.add(KE_current/N);
					U_list.add((PE_current + KE_current)/N);
					temp_list.add(0.66667*KE_current/N);
					data.println(time_current + " " + PE_current + " " + KE_current + " " + (PE_current+KE_current) + " " + (0.66667*KE_current/N));
				}
				PrintWriter traj = null;
				try {
					traj = new PrintWriter(filePath + "MD_traj" + Double.toString(kBT) + ".txt", "UTF-8");
				}
				catch (FileNotFoundException | UnsupportedEncodingException e) {
					e.printStackTrace();
				}
				for (int i = 0; i < ALL_PARTICLES.length; i++) {
					traj.println(ALL_PARTICLES[i].x_coordinate + " " + ALL_PARTICLES[i].y_coordinate + " " + ALL_PARTICLES[i].orientation + " " + ALL_PARTICLES[i].ID + " " + ALL_PARTICLES[i].velocity[0] + " " + ALL_PARTICLES[i].velocity[1] + " " + ALL_PARTICLES[i].velocity[2]);
				}
		        traj.print(q[0] + " " + q[1]);
				traj.close();
			} // check step
			
			time_current += time_step;
			currentStep ++;
		}
		
		// Data Analysis
		double sdfn = 0.0;
		double sdfc = 0.0;
		double sdfs = 0.0;
		if (msdn_list.size() > 1) {
			double sef_time = time_total - time_step*eqStep2;
			double rise = (msdn_list.get(msdn_list.size() - 1) - msdn_list.get(0))/n; // total distance traveled per particle
			sdfn = rise/(4*sef_time);
		}
		if (msdc_list.size() > 1) {
			double sef_time = time_total - time_step*eqStep2;
			double rise = (msdc_list.get(msdc_list.size() - 1) - msdc_list.get(0))/c; // total distance traveled per particle
			sdfc = rise/(4*sef_time);
		}
		if (msds_list.size() > 1) {
			double sef_time = time_total - time_step*eqStep2;
			double rise = (msds_list.get(msds_list.size() - 1) - msds_list.get(0))/s; // total distance traveled per particle
			sdfs = rise/(4*sef_time);
		}
		if (rdf_list.size() > 0) {
			double[] radial_av = MC_utils.list_average(rdf_list); // average rdf over each data collection
			double[] radial_box = MC_utils.boxcar(radial_av, 5); // should be 5*N - 4 long
			PrintWriter radial = null;
			try {
				radial = new PrintWriter(filePath + "MD_radial" + Double.toString(kBT) + ".txt", "UTF-8");
			}
			catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < radial_box.length; i++) {
				double distance = i*Math.min(q[0]/2.0, q[1]/2.0) / (5*N);
				radial.println(distance + " " + radial_box[i]);
			}
			radial.close();
		}
		
		// Data Report
		System.out.println("\n\nData Report:");
		System.out.println("MB Waters: " + n + " (sigma: " + sigma_n + " epsilon: " + epsilon_n + ")");
		System.out.println("LJ_c particles: " + c + " (sigma: " + sigma_c + " epsilon: " + epsilon_c + ")");
		System.out.println("LJ_s particles: " + s + " (sigma: " + sigma_s + " epsilon: " + epsilon_s + ")");
		System.out.println("----------------------------------------------------");
		System.out.println("Set T*: " + kBT + " (tau: " + tau + ")");
		System.out.println("V*: " + Math.round(100*q[0]*q[1])/100.0 + " (rho: " + rho + ")");
		System.out.println("Time: " + time_total + " (dt: " + time_step + ")");
		System.out.println("eq. Time 1: " + Double.toString(eqStep1*time_step));
		System.out.println("eq. Time 2: " + Double.toString(eqStep2*time_step));
		System.out.println("----------------------------------------------------");
		if (U_list.size() > 1) {
			System.out.println("PE: " + MC_utils.stats(PE_list));
			System.out.println("KE: " + MC_utils.stats(KE_list));
			System.out.println("Total energy: " + MC_utils.stats(U_list));
			System.out.println("T*: " + MC_utils.stats(temp_list));
			System.out.println("----------------------------------------------------");
		}
		if (pres_list.size() > 1) {
			System.out.println("P*: " + MC_utils.stats(pres_list));
			System.out.println("Heat Capacity: " + MC_utils.stats(heat_list));
			System.out.println("Self-Diffusion Coefficients: (MB: " + Double.toString(sdfn) + ", LJ_c: " + Double.toString(sdfc) + ", LJ_s: " + Double.toString(sdfs) + ")");
			System.out.println("----------------------------------------------------");
		}
		System.out.println(MC_utils.time_message((System.currentTimeMillis() - start_time)));
		System.out.print("----------------------------------------------------\n");
		
		// Save coordinates to file
		PrintWriter pic = null;
		try {
			pic = new PrintWriter(filePath + "MD_pic" + Double.toString(kBT) + ".txt", "UTF-8");
		}
		catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		for (int i = 0; i < ALL_PARTICLES.length; i++) {
			pic.println(ALL_PARTICLES[i].x_coordinate + " " + ALL_PARTICLES[i].y_coordinate + " " + ALL_PARTICLES[i].orientation + " " + ALL_PARTICLES[i].ID);
		}
        pic.print(q[0] + " " + q[1]);
		pic.close();
		data.close();
	}
}
