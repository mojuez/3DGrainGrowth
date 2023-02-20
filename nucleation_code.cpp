#include<iostream>
#include<cstdlib>
#include<math.h>
#include<string>
#include<iostream>
#include<fstream>
#include<chrono>
#include<random>


int Nx = 128, Ny = 128, Nz = 128; // system dimension
double dx = 1.0, dy = 1.0, dz = 1.0;
int Norient = 36; // number of orientations
double dt = 1e-7;

double pi = 3.14159265358979323846264;
double kappac = 0.1;
//double sigma = sqrt(0.5 * kappac)/3.0/1.2;
double lamda = 1.0;

double Ttrans = 200; //K--Temperature
double Ttemp = 300;
double Q = 7e5;
int distance = 6;

//nucelation rate parameters
double N = 0.01;
double kb = 1.3806452 * pow(10, -23); //boltzmann constant J / K.




//3D index to 1D index
int convert3Dindex(int i, int j, int k) {
	int index1D;
	index1D = i * Ny * Nz + j * Nz + k;
	return index1D;
}

//4D index to 1D index
int convert4Dindex(int i, int j, int k, int n) {
	int index1D;
	index1D = i * Ny * Nz * Norient + j * Nz * Norient + k * Norient + n;
	return index1D;
}

// initialization of eta
void initialize(double* eta) {
	srand(30);
	for (int i = 0; i < Nx * Ny * Nz * Norient; i++) {
		eta[i] = ((double)rand() / RAND_MAX) * 0.002 - 0.001; // eta in the range of (-0.001, 0.001)
	}
}

void output_eta(double* eta, int step) {
	int index_eta;
	std::string filename = "eta" + std::to_string(step) + ".txt";
	std::ofstream outputfile;
	outputfile.open(filename);
	//outputfile << "x" << " " << "y" << " " << "z" << " " << "eta" << std::endl;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				outputfile << i << " " << j << " " << k << " ";
				for (int n = 0; n < Norient; n++) {
					index_eta = convert4Dindex(i, j, k, n);
					outputfile << eta[index_eta] << " ";
				}				
				outputfile << std::endl;
			}
		}
	}
}

double outputgrainvol(double* eta, double* grainvol, int step) {
	int indexeta, indexgrainvol;
	double avegrainvol = 0;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				// get corre
				indexgrainvol = convert3Dindex(i, j, k);
				grainvol[indexgrainvol] = 0;
				for (int n = 0; n < Norient; n++) {
					indexeta = convert4Dindex(i, j, k, n);
					grainvol[indexgrainvol] += pow(eta[indexeta], 2);
				}
				// calcualate the average grainvol
				avegrainvol += grainvol[indexgrainvol];
			}
		}
	}
	avegrainvol = avegrainvol / (Nx * Ny * Nz);
	// output the grainvol
	std::string filename = "grainvol" + std::to_string(step) + ".txt";
	std::ofstream outputfile;
	outputfile.open(filename);
	outputfile << "x" << " " << "y" << " " << "z" << " " << "grainvol" << std::endl;
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				indexgrainvol = convert3Dindex(i, j, k);
				outputfile << i << " " << j << " " << k << " " << grainvol[indexgrainvol] << std::endl;
			}
		}
	}
	outputfile.close();

	return avegrainvol;

}


void get_nucleation_growth(double* eta) {
	int index1D_4D;
	double sum_eta_pow2 = 0; // sum of eta**2
	//double sum_eta_pow2_2 = 0;
	double total_energy = 0; // total energy = dgv
	double volume = Nx * Ny * Nz * dx * dy * dz;
	double random;
	double A, B, C;
	double dgm;

	//nucleation energy related need to change
	if (Ttemp <= Ttrans) {
		A = 0.3 * Q;
	}
	else {
		A = (0.8 + 0.06 * (Ttemp - Ttrans)) * Q;
	}
	dgm = Q * (Ttemp - Ttrans) / Ttrans;
	B = 3 * A - 12 * dgm;
	C = 2 * A - 12 * dgm;

	for (int i = 0; i < Nx * Ny * Nz * Norient; i++) {
		//index1D_4D = convert4Dindex(i, j, k, n);
		//sum_eta_pow2 += pow(eta[index1D_4D], 2);
		//total_energy += A / 2 * pow(eta[index1D_4D], 2) - B / 3 * pow(eta[index1D_4D], 3);
		sum_eta_pow2 += pow(eta[i], 2);
		total_energy += A / 2 * pow(eta[i], 2) - B / 3 * pow(eta[i], 3);
	}
	//total_energy += C / 4 * pow(sum_eta_pow2, 2);
	total_energy += C / 4 * pow(sum_eta_pow2, 2);

	//std::cout << total_energy << std::endl;

	/*double dgv = total_energy;
	double dgs = 16 * pi * pow(2 * kappac / (3 * lamda), 3) / (3 * pow(dgv, 2));
	double nucleation_rate = N * exp(-dgs / (kb * Ttemp));
	double nucleation_probability = 1.0 - exp(-nucleation_rate * dt);
	
	double radius_critic = -4.0 * kappac / (3 * lamda * dgv);;*/

	double nucleation_probability = 0.4;
	double radius_critic = 2;

	//std::cout << nucleation_rate << std::endl;
	std::cout << nucleation_probability << std::endl;
	std::cout << radius_critic << std::endl;

	srand(50);
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {
			for (int k = 0; k < Nz; k++) {
				//index1D_3D = convert3Dindex(i, j, k);


				random = ((double)rand() / RAND_MAX); //range 0~1

				sum_eta_pow2 = 0;

				for (int n = 0; n < Norient; n++) {
					index1D_4D = convert4Dindex(i, j, k, n);
					sum_eta_pow2 += pow(eta[index1D_4D], 2);
				}

				//change to sum eta**2 as judgement
				if (random < nucleation_probability && sum_eta_pow2 < 0.9) {
					//calculate the critical radius
					

					if ((i - radius_critic - distance) >= 1 && (i + radius_critic + distance) <= Nx &&
						(j - radius_critic - distance) >= 1 && (j + radius_critic + distance) <= Ny &&
						(k - radius_critic - distance) >= 1 && (k + radius_critic + distance) <= Nz) {
						int overlap = 0;

						for (int ii = i - radius_critic - distance; ii < i + radius_critic + distance; ii++) {
							for (int jj = j - radius_critic - distance; jj < j + radius_critic + distance; jj++) {
								for (int kk = k - radius_critic - distance; kk < k + radius_critic + distance; kk++) {
									if (pow((ii - i), 2) + pow((jj - j), 2) + pow((kk - k), 2)
										<= pow(radius_critic + distance, 2)) {

										sum_eta_pow2 = 0;
										//add sum eta**2 > 0.9 for loop grain overlap
										for (int nn = 0; nn < Norient; nn++) {
											index1D_4D = convert4Dindex(ii, jj, kk, nn);
											sum_eta_pow2 += pow(eta[index1D_4D], 2);
										}
										if (sum_eta_pow2 > 0.9) {
											overlap = overlap + 1;
										}
									}
								}
							}
						}

						if (overlap == 0) {
							int random_2 = 35 * ((double)rand() / RAND_MAX); //range 0~35 * numer pf orientation = 36
							for (int ii = i - radius_critic - distance; ii < i + radius_critic + distance; ii++) {
								for (int jj = j - radius_critic - distance; jj < j + radius_critic + distance; jj++) {
									for (int kk = k - radius_critic - distance; kk < k + radius_critic + distance; kk++) {
										if (pow((ii - i), 2) + pow((jj - j), 2) + pow((kk - k), 2)
											<= pow(radius_critic, 2)) {

											//for loop  eta[random_2] = 1, others = 0
											for (int nn = 0; nn < Norient; nn++) {
												index1D_4D = convert4Dindex(ii, jj, kk, nn);
												if (nn == random_2) {
													eta[index1D_4D] = 1;
												}
												else {
													eta[index1D_4D] = 0;
												}
											}
										}
									}
								}
							}
							std::cout << i << "" << j << "" << k << std::endl;
						}
					}
				}
			}
		}
	}
}






int main() {
	double* eta = new double[Nx * Ny * Nz * Norient]; // eta
	double* grainvol = new double[Nx * Ny * Nz];
	//assign random number to eta
	initialize(eta);

	outputgrainvol(eta, grainvol, 0);
	
	//read eta and inset into the funtion get_nucleation_growth
	get_nucleation_growth(eta);
	output_eta(eta, 1);
	outputgrainvol(eta, grainvol, 1);

}

