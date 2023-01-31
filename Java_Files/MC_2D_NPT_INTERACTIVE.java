package mercedes_benz;

/* Daniel Gray
 * Goal: copy 2D MC code I have in python to make it faster
 * 
 * Change Log:
 * 16-June-2018:
 * Created file
 * 6-July-2018:
 * Made it working
 * Tested and it converges to same energy
 * 11-July-2018:
 * Added msd, rdf, and saving all data to files
 * 12-July-2018:
 * Finished rdf, added pressure
 * 17-July-2018:
 * Added density quadratic
 * Made box rectangular by default
 * 23-July-2018:
 * Removed density quadratic because it cause issues
 * Fixed pressure issue
 * 
 * Pseudocode:
 * Parameters
 * Initialization
 * Main Loop
 * Data Analysis
 * Data Report
 */

import java.util.concurrent.ThreadLocalRandom;
//import java.util.Arrays;
import java.util.ArrayList;
import java.util.Arrays;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

public class MC_2D_NPT {

	public static void main(String[] args) {
		
		int n = 0; // MB waters
		int c = 0; // small crowders
		int s = 0; // large crowders
		double kBT = 0.0;
		double pressure = 0.0;
		int PASSES = 0; // these are determined when program ran
		boolean restart = false;
		n = MC_utils.choose(n, "Enter the number of MB waters:");
		c = MC_utils.choose(c, "Enter the number of large LJ particles:");	
		s = MC_utils.choose(s, "Enter the number of small LJ particles:");	
		kBT = MC_utils.choose(kBT, "Enter the temperature:");
		pressure = MC_utils.choose(pressure, "Enter the pressure: ");
		PASSES = MC_utils.choose(PASSES, "Enter the number of passes:");
		restart = MC_utils.choose(restart, "Press 'y' for restart (anything else to not):");
		
		final long start_time = System.currentTimeMillis();	
		
		// Parameters
		final int N = n + c + s; // total particles
		int STEPS = PASSES*N; // steps to run
		double density = 1.0; // initial density of water
		
		final String filePath = "C:\\Users\\Yo-Dan\\Documents\\School\\Internships\\Summer 2018-Duquesne\\Java_MC_Data\\";
		
		final double sigma_c = 50; // LJ parameters for crowding particles
		final double epsilon_c = 0.1;
		final double sigma_s = 0.7;
		final double epsilon_s = 0.1;
		
		final double Lx = Math.sqrt(N/density); // rectangle
		final double Ly = Math.sqrt(N/density);
		double[] q = {Lx, Ly};
		final double sigma_MB = 0.7; // LJ parameters for MB water
		final double epsilon_MB = 0.1;
		final double r_HB = 1.0; // HB parameters for MB water
		final double sigma_HB = 0.085;
		final double epsilon_HB = -1.0;
		
		int eqPass1 = 250000; // energy/density
		int eqPass2 = 9000000; // everything else
		if (PASSES < eqPass1) {
			eqPass1 = 0;
			eqPass2 = 0;
		}
		else if (PASSES < eqPass2) {
			eqPass2 = eqPass1;
		}
		final int dataFreq = 1000*N; // data saved into lists
		final int printFreq = 200000*N; // how often status printed to screen
		final int volumeFreq = 5*N; // volume moves
		final int ddqFreq = 100*volumeFreq; // dq changes
		final int ddsFreq = 400; // ds changes
		final int ddoFreq = (int)(ddsFreq*1.5); // do changes
		
		final double[] sigma = new double[N]; // LJ terms
		final double[] epsilon = new double[N];
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
				sigma[i] = sigma_MB;
				epsilon[i] = epsilon_MB;
			}
		}
		final double G_constant = -1.0/(2*sigma_HB*sigma_HB); // For Gaussian
		ArrayList<Double> energy_list = new ArrayList<> (); // data lists
		ArrayList<Double> density_list = new ArrayList<> ();
		ArrayList<Double> heat_list = new ArrayList<> ();
		ArrayList<Double> isothermal_list = new ArrayList<> ();
		ArrayList<Double> expansion_list = new ArrayList<> ();
		ArrayList<double[]> radial_list = new ArrayList<> ();
		ArrayList<Double> msd_list = new ArrayList<> ();
		double msd = 0;
		
		final double accRatio = 0.5; // acceptance ratio
		double delq = 0.01*Math.min(q[0], q[1]); // max moves
		double delo = 0.5;
		double dels = 0.1;
		final double ddq = 0.05*delq; // max change in the max moves
		final double ddo = 0.05*delo;
		final double dds = 0.05*dels;
		double qAccRatio = 0; // running acceptance frequencies
		double oAccRatio = 0;
		double sAccRatio = 0;
		
		int qMove = 0; // counters
		int oMove = 0;
		int sMove = 0;
		int qAccept = 0;
		int oAccept = 0;
		int sAccept = 0;
		int currentStep = 1;
		
		// Program execution
		Particle[] ALL_PARTICLES = new Particle[N];
		if (restart) {
			String filename = "";
			filename = MC_utils.choose("Enter name of structure file");
			System.out.println("Writing in particles.");
			ALL_PARTICLES = MC_utils.restart(n, c, s, filePath + filename + ".txt", q);
		}
		else {
			System.out.println("Initializing particles.");
			ALL_PARTICLES = MC_utils.setup_particles(n, s, c, q, sigma);
		}
		
		double Ui = MC_utils.U_total(ALL_PARTICLES, q, sigma, epsilon, G_constant, r_HB, epsilon_HB);
//		MC_utils.print_particles(ALL_PARTICLES);
		System.out.println("Equilibrating. Ui = " + Ui);
		System.out.println(Arrays.toString(q));
		
		PrintWriter data = null;
		try {
			data = new PrintWriter(filePath + "MC_data" + Double.toString(kBT) + ".txt", "UTF-8");
		}
		catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}
		
		while (currentStep <= STEPS) {
			if (currentStep%volumeFreq == 0) { // volume move
				double vOld = q[0]*q[1];
				double amount = (2*Math.random() - 1)*delq;
				double x_scale = (amount + q[0])/q[0];
				double y_scale = (amount + q[1])/q[1];
				q[0] *= x_scale; // change box size
				q[1] *= y_scale;
				for (int i = 0; i < N; i++) { // change particle positions
					ALL_PARTICLES[i].x_coordinate *= x_scale;
					ALL_PARTICLES[i].y_coordinate *= y_scale;
				}
				double Uf = MC_utils.U_total(ALL_PARTICLES, q, sigma, epsilon, G_constant, r_HB, epsilon_HB);
				double weight = (Uf-Ui) + pressure*(q[0]*q[1]-vOld) - kBT*N*Math.log(q[0]*q[1]/vOld); // dH = dU + PdV - kTNln(V'/V)
				if (MC_utils.Accept(weight, kBT)) {
					Ui = Uf;
					qAccept ++;
				}
				else {
					q[0] /= x_scale; // reset
					q[1] /= y_scale;
					for (int i = 0; i < N; i++) {
						ALL_PARTICLES[i].x_coordinate /= x_scale;
						ALL_PARTICLES[i].y_coordinate /= y_scale;
					}
				}
				qMove ++;
			}
			else {
				int particle = ThreadLocalRandom.current().nextInt(0, N); // random particle
				int parameter;
				if (particle < n) {
					parameter = ThreadLocalRandom.current().nextInt(0, 3); // random parameter
				}
				else {
					parameter = ThreadLocalRandom.current().nextInt(0, 2); // no orientation for LJ particles
				}
				if (parameter == 2) { // orientation move
					double amount = (2*Math.random() - 1)*delo;
					double U_i = MC_utils.U_particle(ALL_PARTICLES, particle, q, sigma, epsilon, G_constant, r_HB, epsilon_HB);
					ALL_PARTICLES[particle].orientation += amount;
					double U_f = MC_utils.U_particle(ALL_PARTICLES, particle, q, sigma, epsilon, G_constant, r_HB, epsilon_HB);
					if (MC_utils.Accept(U_f-U_i, kBT)) {
						Ui += U_f - U_i;
						oAccept ++;
					}
					else { // move rejected
						ALL_PARTICLES[particle].orientation -= amount;
					}
					oMove ++;
				}
				else { // translation move
					double amount = (2*Math.random()-1)*dels;
					double U_i = MC_utils.U_particle(ALL_PARTICLES, particle, q, sigma, epsilon, G_constant, r_HB, epsilon_HB);
					if (parameter == 0) {
						ALL_PARTICLES[particle].x_coordinate += amount;
					}
					else {
						ALL_PARTICLES[particle].y_coordinate += amount;
					}
					double U_f = MC_utils.U_particle(ALL_PARTICLES, particle, q, sigma, epsilon, G_constant, r_HB, epsilon_HB);
					if (MC_utils.Accept(U_f-U_i, kBT)) {
						Ui += U_f - U_i;
						msd += amount*amount; // add to msd variable
						if (parameter == 0) { // wrap particle here
							if (ALL_PARTICLES[particle].x_coordinate < 0) {								
								ALL_PARTICLES[particle].x_coordinate += q[0];
							}
							else if (ALL_PARTICLES[particle].x_coordinate > q[0]) {
								ALL_PARTICLES[particle].x_coordinate -= q[0];
							}
						}
						else { // wrap particle here
							if (ALL_PARTICLES[particle].y_coordinate < 0) {								
								ALL_PARTICLES[particle].y_coordinate += q[1];
							}
							else if (ALL_PARTICLES[particle].y_coordinate > q[1]) {
								ALL_PARTICLES[particle].y_coordinate -= q[1];
							}						}
						sAccept ++;
					}
					else { // move rejected
						if (parameter == 0)
							ALL_PARTICLES[particle].x_coordinate -= amount;
						else
							ALL_PARTICLES[particle].y_coordinate -= amount;
					}
					sMove ++;
				}
			}
			
			if (currentStep%ddqFreq == 0) { // Change volume move size
				qAccRatio = (float)qAccept/(float)qMove;
				if (qAccRatio > (accRatio + 0.05)) {
					delq += ddq;
				}
				if ((qAccRatio < (accRatio - 0.05)) && (delq > (3*ddq/2))) {
					delq -= ddq;
				}
				qAccept = 0;
				qMove = 0;
			}
			if (currentStep%ddoFreq == 0) { // Change orientation move size
				oAccRatio = (float)oAccept/(float)oMove;
				if (oAccRatio > (accRatio + 0.05)) {
					delo += ddo;
				}
				if ((oAccRatio < (accRatio - 0.05)) && (delo > (3*ddo/2))) {
					delo -= ddo;
				}
				oAccept = 0;
				oMove = 0;
			}
			if (currentStep%ddsFreq == 0) { // Change translation move size
				sAccRatio = (float)sAccept/(float)sMove;
				if (sAccRatio > (accRatio + 0.05)) {
					dels += dds;
				}
				if ((sAccRatio < (accRatio - 0.05)) && (dels > (3*dds/2))) {
					dels -= dds;
				}
				sAccept = 0;
				sMove = 0;
			}
			
			if (currentStep%dataFreq == 0) { // Save data to lists
				Ui = MC_utils.U_total(ALL_PARTICLES, q, sigma, epsilon, G_constant, r_HB, epsilon_HB); // to avoid long-term roundoff errors
				if (currentStep%printFreq == 0) {
					System.out.println("Pass " + (int)(currentStep/N));
				}
				if (currentStep > eqPass2*N) {
					energy_list.add(Ui/N);
					density_list.add(N/(q[0]*q[1])); // after equilibrating just keep extending list
					double h = MC_utils.C_p(energy_list, N, kBT);
					double i = MC_utils.kap(density_list, N, kBT); // equilibrate thermodynamics
					double e = MC_utils.alpha(energy_list, density_list, N, kBT);
					heat_list.add(h);
					isothermal_list.add(i);
					expansion_list.add(e);
					radial_list.add(MC_utils.radDist(ALL_PARTICLES, 5*N, Math.min(q[0]/2.0, q[1]/2.0), q));
					msd_list.add(msd/N);
					data.println(currentStep/N + " " + Ui/N + " " + N/(q[0]*q[1]) + " " + h + " " + i + " " + e + " " + msd/N);
				}
				else if (currentStep > eqPass1*N) {
					energy_list.add(Ui/N);
					density_list.add(N/(q[0]*q[1])); // after equilibrating just keep extending list
					data.println(currentStep/N + " " + Ui/N + " " + N/(q[0]*q[1]));
				}
				PrintWriter traj = null;
				try {
					traj = new PrintWriter(filePath + "MC_traj" + Double.toString(kBT) + ".txt", "UTF-8");
				}
				catch (FileNotFoundException | UnsupportedEncodingException e) {
					e.printStackTrace();
				}
				for (int i = 0; i < N; i++) {
					traj.println(ALL_PARTICLES[i].x_coordinate + " " + ALL_PARTICLES[i].y_coordinate + " " + ALL_PARTICLES[i].orientation + " " + ALL_PARTICLES[i].ID);
				}
		        traj.print(q[0] + " " + q[1]);
				traj.close();
			}
			
			currentStep++;
		} // main loop
		
		// data analysis
		double sdf = 0;
		if (msd_list.size() > 1) {
			double msd_passes = (PASSES - eqPass2)*0.666667*(volumeFreq-1)/volumeFreq; // correct number of passes (run)
			double rise = (msd_list.get(msd_list.size() - 1) - msd_list.get(0))/N; // total distance traveled per particle
			sdf = rise/(4*msd_passes); // self-diffusion coefficient = msd_slope/(2*d) with d = dimensionality = 2
		}
		if (radial_list.size() > 0) {
			double[] radial_av = MC_utils.list_average(radial_list); // average rdf over each data collection
			double [] radial_box = MC_utils.boxcar(radial_av, 5); // should be 5*N - 4 long
			PrintWriter radial = null;
			try {
				radial = new PrintWriter(filePath + "MC_radial" + Double.toString(kBT) + ".txt", "UTF-8");
			}
			catch (FileNotFoundException | UnsupportedEncodingException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < radial_box.length; i++) {
				double distance = i*Math.min(q[0]/2.0, q[1]/2.0 / (5*N));
				radial.println(distance + " " + radial_box[i]); // just a graph
			}
			radial.close();
		}
		
		// Data Report
		System.out.println("\n\nData Report:");
		System.out.println("MB Waters: " + n + " (sigma: " + sigma_MB + " epsilon: " + epsilon_MB + ")");
		System.out.println("LJ_C particles: " + c + " (sigma: " + sigma_c + " epsilon: " + epsilon_c + ")");
		System.out.println("LJ_S particles: " + s + " (sigma: " + sigma_s + " epsilon: " + epsilon_s + ")");
		System.out.println("T* = " + kBT);
		System.out.println("P* = " + pressure);
		System.out.println("Steps: " + STEPS + " (Passes: " + STEPS/N + ")");
		System.out.println("eq. Pass 1: " + Integer.toString(eqPass1));
		System.out.println("eq. Pass 2: " + Integer.toString(eqPass2));
		System.out.println("-----------------------------------------------");
		if (energy_list.size() > 1) {
			System.out.println("U* = " + MC_utils.stats(energy_list));
			System.out.println("Density = " + MC_utils.stats(density_list));
			System.out.println("-----------------------------------------------");
		}
		System.out.println("Volume AccRatio = " + Math.round(qAccRatio*1000)/1000.0 + " (delq = " + delq + ")");
		System.out.println("Orientation AccRatio = " + Math.round(oAccRatio*1000)/1000.0 + " (delo = " + delo + ")");
		System.out.println("Translation AccRatio = " + Math.round(sAccRatio*1000)/1000.0 + " (dels = " + dels + ")");
		System.out.println("-----------------------------------------------");
		if (heat_list.size() > 1) {
			System.out.println("heat capacity = " + MC_utils.stats(heat_list));
			System.out.println("isothermal compressibility = " + MC_utils.stats(isothermal_list));
			System.out.println("thermal expansion coefficient = " + MC_utils.stats(expansion_list));
			System.out.println("self-diffusion coefficient = " + Double.toString(sdf));
			System.out.println("-----------------------------------------------");
		}
		System.out.println(MC_utils.time_message((System.currentTimeMillis() - start_time)));
		System.out.print("-----------------------------------------------\n");
		
		// Save coordinates to file
		PrintWriter pic = null;
		try {
			pic = new PrintWriter(filePath + "MC_pic" + Double.toString(kBT) + ".txt", "UTF-8");
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
