package mercedes_benz;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class MC_utils {
	
	public static int choose(int number, String name) { // allows user input for a int variable
		Scanner input = new Scanner(System.in);
		System.out.println(name);
		while (input.hasNextLine()) {
			String num = input.nextLine();
			try {
				number = Integer.parseInt(num);
				break;
			}
			catch (NumberFormatException e) {
				System.out.println("Incorrect input.  Enter an integer.");
			}
		}
		return number;
	}
	public static double choose(double number, String name) {  // choose a double variable
		Scanner input = new Scanner(System.in);
		System.out.println(name);
		while (true) {
			String num = input.nextLine();
			try {
				number = Double.parseDouble(num);
				break;
			}
			catch (NumberFormatException e) {
				System.out.println("Incorrect input.  Enter a float.");
			}
		}
		return number;
	}
	public static boolean choose(boolean s, String name) { // choose a boolean
		Scanner input = new Scanner(System.in);
		System.out.println(name);
		String answer = input.nextLine();
		if (answer.equalsIgnoreCase("y")) {
			s = true;
		}
		return s;
	}
	public static String choose(String name) { // choose a string
		Scanner input = new Scanner(System.in);
		System.out.println(name);
		String answer = input.nextLine();
		return answer;
	}
	public static Particle[] setup_particles(int n, int s, int c, double[] q, double[] sigma) { // randomly sets up particles without allowing overlap
		Particle[] AP = new Particle[(n+s+c)];
		for (int i = 0; i < (n+s+c); i++) {
			Particle particle = new Particle();
			if (i < c) {
				particle.setID("LJ_c"); // large LJ particles first
				particle.setO(0.0);
			}
			else if (i < (c+s)) { // small LJ particle
				particle.setID("LJ_s");
				particle.setO(0.0);
			}
			else {
				particle.setID("MB"); // water
				particle.setO(Math.random()*Math.PI);
			}
			boolean tooClose = true;
			int tries = 0;
			double[] bestTry = {0,0,0};
			while (tooClose) {
				tooClose = false;
				particle.setX(Math.random()*q[0]); // put particle in box randomly
				particle.setY(Math.random()*q[1]);
				for (int j = 0; j < i; j++) { // test if overlaps with any other particle
					double[] u_hat_not = nearest_image(particle, AP[j], q);
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
				if (tries >= 10000*(n+s+c)) { // give up after 10000 passes
					particle.setX(bestTry[0]);
					particle.setY(bestTry[1]);
					System.out.println("Particle too close.");
					break;
				}
			}
			AP[i] = particle; // append particle
		}
		return AP;
	}
	public static void print_particles(Particle[] AP) { // outputs particle information
		for (int i = 0; i < AP.length; i++) {
			String[] particle = new String[4];
			particle[0] = AP[i].getID();
			particle[1] = Double.toString(Math.round(AP[i].getX()*1000)/1000.0);
			particle[2] = Double.toString(Math.round(AP[i].getY()*1000)/1000.0);
			particle[3] = Double.toString(Math.round(AP[i].getO()*1000)/1000.0);
			System.out.println(Arrays.toString(particle));
		}
	}
	public static Particle[] restart(int n, int c, int s, String name, double[] q) { // start with input structure
		Particle[] AP = new Particle[n+c+s];
		FileReader position = null;
		try {
			position = new FileReader(name);
			BufferedReader pos = new BufferedReader(position);
			for (int i = 0; i < (n+c+s); i++) {
				Particle particle = new Particle();
				String line = pos.readLine();
				String[] ln = line.split(" ");
				double x = 0;
				double y = 0;
				double o = 0;
				String id = "";
				if (ln.length >= 3) {
					try {
						x = Double.parseDouble(ln[0]);
						y = Double.parseDouble(ln[1]);
						o = Double.parseDouble(ln[2]);
					}
					catch (NumberFormatException e) {
						System.out.println("Error on particle " + Integer.toString(i)); // if wrong type of file
					}
				}
				if (i < c) {
					id = "LJ_c";
				}
				else if (i < (c+s)) {
					id = "LJ_s";
				}
				else {
					id = "MB";
				}
				particle.setX(x);
				particle.setY(y);
				particle.setO(o);
				particle.setID(id);
				AP[i] = particle;
			}
			String line = pos.readLine();
			String[] ln = line.split(" ");
			q[0] = Double.parseDouble(ln[0]);
			q[1] = Double.parseDouble(ln[1]);
			pos.close();
			position.close();
		}
		catch (IOException e) {
			System.out.println("Couldn't open or read file."); // if incorrect file name
			e.printStackTrace();
		}
		return AP;
	}
	
	public static double[] nearest_image(Particle i, Particle j, double[] q) { // nearest image across periodic box
		double dx = j.getX() - i.getX();
		double dy = j.getY() - i.getY(); // initial distances apart
		if (dx > (q[0]/2)) { // test if closer across box boundaries
			dx -= q[0];
		}
		if (dx < (-q[0]/2)) {
			dx += q[0];
		}
		if (dy > (q[1]/2)) {
			dy -= q[1];
		}
		if (dy < (-q[1]/2)) {
			dy += q[1];
		}
		double[] u_hat_not = {dx, dy};
		return u_hat_not;
	}
	
	public static double U_total(Particle[] AP, double[] q, double[] sigma, double[] epsilon, double G, double rHB, double epHB) {
		double U = 0;
		for (int i = 0; i < AP.length; i++) {
			for (int j = (i+1); j < AP.length; j++) { // find each pair
				double sigma_av = (sigma[i] + sigma[j]) / 2;
				double epsilon_av = (epsilon[i] + epsilon[j]) / 2;
				U += U_pair(AP[i], AP[j], q, sigma_av, epsilon_av, G, rHB, epHB);
			}
		}
		return U;
	}
	public static double U_pair(Particle i, Particle j, double[] q, double sigma, double epsilon, double G, double rHB, double epHB) {
		double U = 0;
		double[] u_hat_not = nearest_image(i, j, q);
		double r_ij = Math.sqrt(u_hat_not[0]*u_hat_not[0] + u_hat_not[1]*u_hat_not[1]);
		if (r_ij < 4*sigma) {
			U += U_LJ(sigma, epsilon, r_ij); // all particles have LJ
		}
		if ((i.getID() == "MB") && (j.getID() == "MB") && (r_ij < 1.5*rHB) && (r_ij > 0.5*rHB)) { // only MB waters H-bond
			double[] u_hat = {u_hat_not[0]/r_ij, u_hat_not[1]/r_ij};
			double[][] i_arms = get_arms(i.getO());
			double[][] j_arms = get_arms(j.getO());
			double i_dot = Math.max(DOT(i_arms[0], u_hat), Math.max(DOT(i_arms[1], u_hat), DOT(i_arms[2], u_hat)));
			double j_dot = Math.min(DOT(j_arms[0], u_hat), Math.min(DOT(j_arms[1], u_hat), DOT(j_arms[2], u_hat)));
			U += U_HB(G, rHB, epHB, u_hat, i_dot, j_dot, r_ij);
		}
		return U;
	}
	public static double U_particle(Particle[] AP, int i, double[] q, double[] sigma, double[] epsilon, double G, double rHB, double epHB) { // energy of one particle
		double U = 0;
		for (int j = 0; j < AP.length; j++) { // find each other particle
			if (i != j) {
				double sig_av = (sigma[i]+sigma[j])/2.0;
				double ep_av = (epsilon[i]+epsilon[j])/2.0;
				U += U_pair(AP[i], AP[j], q, sig_av, ep_av, G, rHB, epHB);
			}
		}
		return U;
	}
	public static double U_LJ(double sigma, double epsilon, double r) {
		double U1 = sigma/r;
		double U2 = U1*U1*U1;
		double U3 = U2*U2;
		double U = 4*epsilon*U3*(U3 - 1); // U = 4*e*((sigma/r)^12 - (sigma/r)^6)
		return U;
	}
	public static double U_HB(double G_constant, double rHB, double epHB, double[] u_hat, double i, double j, double r) {
		return epHB*Gaussian(rHB - r, G_constant)*Gaussian(i - 1, G_constant)*Gaussian(j + 1, G_constant); // i and j are max and min dot products
	}
	public static double[][] get_arms(double angle) { // returns vector of arm vectors
		double[] arm1 = {Math.cos(angle), Math.sin(angle)};
		double[] arm2 = {Math.cos(angle + (2*Math.PI/3.0)), Math.sin(angle + (2*Math.PI/3.0))};
		double[] arm3 = {Math.cos(angle - (2*Math.PI/3.0)), Math.sin(angle - (2*Math.PI/3.0))};
		double[][] arms = {arm1, arm2, arm3};
		return arms;
	}
	public static double Gaussian(double x, double G) {
		return Math.exp(x*x*G);
	}
	
	public static double[] vector_add(double[] list1, double[] list2) { // basic vector stuff
		double[] list3 = new double[list1.length];
		if (list1.length == list2.length) {
			for (int i = 0; i < list1.length; i++) {
				list3[i] = list1[i] + list2[i];
			}
		}
		return list3; // returns blank list if not same size
	}
	public static double[] vector_sub(double[] list1, double[] list2) {
		double[] list3 = new double[list1.length];
		if (list1.length == list2.length) {
			for (int i = 0; i < list1.length; i++) {
				list3[i] = list1[i] - list2[i];
			}
		}
		return list3; // returns blank list if not same size
	}
	public static double[] vector_mult(double[] list, double k) {
		double[] listk = new double[list.length];
		for (int i = 0; i < list.length; i++) {
			listk[i] = list[i]*k;
		}
		return listk;
	}
	public static double DOT(double[] i, double[] j) {
		double result = 0;
		for (int k = 0; k < i.length; k++) {
			result += (i[k]*j[k]);
		}
		return result;
	}
	
	public static boolean Accept(double delH, double kBT) {
		if (delH <= 0) {
			return true; // always accept unless increases energy
		}
		else {
			double probAccept = Math.exp(-delH/kBT); // probability move accepted
			double randDraw = Math.random();
			if (randDraw < probAccept) {
				return true;
			}
			else {
				return false;
			}
		}
	}
	
	public static double C_p(ArrayList<Double> energy, int N, double T) { // heat capacity of system
		double energy_av = 0;
		double squared_av = 0;
		for (double i : energy) {
			energy_av += (N*i);
			squared_av += (N*N*i*i);
		}
		energy_av /= energy.size();
		squared_av /= energy.size();
		return ((squared_av - (energy_av*energy_av))/(N*T*T)); // Cp = ( <H^2> - <H>^2 ) / ( N*T*T)
	}
	public static double kap(ArrayList<Double> density, double N, double T) { // isothermal compressibility
		double volume_av = 0;
		double squared_av = 0;
		for (double i : density) {
			volume_av += i*N;
			squared_av += (i*i*N*N);
		}
		volume_av /= density.size();
		squared_av /= density.size();
		return ((squared_av - (volume_av*volume_av) )/ (T*volume_av)); // k = ( <V^2> - <V>^2 ) / (T*<V>)
	}
	public static double alpha(ArrayList<Double> energy, ArrayList<Double> density, double N, double T) {
		double volume_av = 0;
		double energy_av = 0;
		double mixed_av = 0;
		for (int i = 0; i < energy.size(); i++) {
			volume_av += N/density.get(i);
			energy_av += N*energy.get(i);
			mixed_av += N*N*energy.get(i)/density.get(i);
		}
		volume_av /= energy.size();
		energy_av /= energy.size();
		mixed_av /= energy.size();
		return ((mixed_av - (volume_av*energy_av)) / (T*T*volume_av)); // a = ( <HV> - <H>*<V> ) / (T*T*<V>)
	}
	
	public static double[] radDist(Particle[] AP, int precision, double maxD, double[] q) {
		double[] radDist = new double[precision]; // how many bins
		for (int i = 0; i < AP.length; i++) { // find radial distribution for each particle and average
			for (int j = 0; j < AP.length; j++) {
				if (i != j) {
					double[] u = nearest_image(AP[i], AP[j], q);
					double r = Math.sqrt(u[0]*u[0] + u[1]*u[1]); // distance between particles
					if (r <= maxD) {
                        int index = (int)(Math.floor(r*precision/maxD)); // which bin the distance goes in
					    double area = 6.283185*r*maxD/precision; // maxD/precision is length of 1 bin
					    radDist[index] += 1/area; // make density list
					}
				}
			}
		}
		radDist = vector_mult(radDist, 1.0/AP.length); // average over all particles
		radDist = vector_mult(radDist, q[0]*q[1]/AP.length); // average over density
		return radDist;
	}
	public static double[] list_average(double[][] list) { // for array
		double[] list_av = new double[list[0].length];
		for (int i = 0; i < list[0].length; i++) { // loop through indices
			double sum = 0;
			for (int j = 0; j < list.length; j++) { // loop through each sublist
				sum += list[j][i];
			}
			sum /= list.length;
			list_av[i] = sum;
		}
		return list_av;
	}
	public static double[] list_average(ArrayList<double[]> list) { // for arraylist
		double[] list_av = new double[list.get(0).length];
		for (int i = 0; i < list.get(0).length; i++) { // loop through indices
			double sum = 0;
			for (int j = 0; j < list.size(); j++) { // loop through each sublist
				sum += list.get(j)[i];
			}
			sum /= list.size();
			list_av[i] = sum;
		}
		return list_av;
	}
	public static double[] boxcar(double[] list, int num) { // smooths graph
		double[] box_list = new double[list.length-num+1];
		for (int i = 0; i < (list.length-num+1); i++) {
			double box_sum = 0;
			for (int j = 0; j < num; j++) {
				box_sum += list[i+j];
			}
			box_sum /= num;
			box_list[i] = box_sum;
		}
		return box_list;
	}
	
	public static double average(double[] list) { // array
		double av = 0;
		for (double i : list) {
			av += i;
		}
		av /= list.length;
		return av;
	}
	public static double average(ArrayList<Double> list) { // arraylist
		double av = 0;
		for (double i : list) {
			av += i;
		}
		av /= list.size();
		return av;
	}
	public static double stdev(double[] list) { // array
		double sd = 0;
		double av = average(list);
		for (double i : list) {
			sd += (i-av)*(i-av);
		}
		sd /= (list.length - 1);
		return Math.sqrt(sd);
	}
	public static double stdev(ArrayList<Double> list) { // arraylist
		double sd = 0;
		double av = average(list);
		for (double i : list) {
			sd += (i-av)*(i-av);
		}
		sd /= (list.size() - 1);
		return Math.sqrt(sd);
	}
	public static String stats(double[] list) { // prints message
		String av = Double.toString(average(list));
		String sd = Double.toString(stdev(list));
		String message = av + " (std. dev. = " + sd + ")";
		return message;
	}
	public static String stats(ArrayList<Double> list) {
		String av = Double.toString(average(list));
		String sd = Double.toString(stdev(list));
		String message = av + " (std. dev. " + sd + ")";
		return message;
	}
	
	public static String time_message(long ms) { // prints elapsed time
		String message = "";
		int seconds = (int)ms/1000;
		int minutes = 0;
		int hours = 0;
		int days = 0;
		while (seconds >= 60) {
			minutes ++;
			seconds -= 60;
		}
		while (minutes >= 60) {
			hours ++;
			minutes -= 60;
		}
		while (hours >= 24) {
			days ++;
			hours -= 24;
		}
		if (days > 0) {
			message = "Elapsed time: " + Integer.toString(days) + " days, " + Integer.toString(hours) + " hours";
		}
		else if (hours > 0) {
			message = "Elapsed time: " + Integer.toString(hours) + " hours, " + Integer.toString(minutes) + " minutes";
		}
		else if (minutes > 0) {
			message = "Elapsed time: " + Integer.toString(minutes) + " minutes, " + Integer.toString(seconds) + " seconds";
		}
		else {
			message = "Elapsed time: " + Integer.toString(seconds) + " seconds";
		}
		return message;
	}

}
