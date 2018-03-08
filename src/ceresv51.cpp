// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <vector>
#include <random>
#include <typeinfo>
#include <algorithm>
#include <iostream>
#include <array>        // std::array
#include <chrono> 
#include <math.h>
#include <set>
#include <ctime>
#include <fstream>

using namespace Rcpp;

double sign_func(double x){
  if(x > 0){
    return 1;
  } else if(x < 0){
  	return -1;
  } else {
  	return 0;
  }
}

double cn_func(double x){
  if(x == 0){
    return 1;
  } else {
  	return 0;
  }
}

double un_nan(double x){
	if(std::isnan(x)){
		return(0);
	} else {
		return(x);
	}
}

int is_less_than_zero(Eigen::VectorXd & x){

	for(int i=0; i < x.size(); ++i){
		if(x(i) >= 0){
			return 0;
		}
	}
	return 1;
}

int is_between_one_and_zero(Eigen::VectorXd & x){

	for(int i=0; i < x.size(); ++i){
		if(x(i) >= 1 || x(i) <= 0){
			return 0;
		}
	}
	return 1;
}

double backtrack(Eigen::SparseMatrix<double> & AtA, Eigen::VectorXd & Atb, Eigen::VectorXd & x, Eigen::VectorXd & delta_x, Eigen::VectorXd & grad_f, Eigen::MatrixXd & SmoothMat,
					double alpha, double beta, double f_x, int NSEGMENTS, int NLINES, int NGENES, int NOBS, int CL_PROB_SIZE, double barrier_t, double LAMBDA_G, double LAMBDA_Smooth){
	//
	double t = 1;
	int n_steps = 1;
	//
	double LAMBDA_STABLE = 1e-8;
	if(NLINES == 1){
		LAMBDA_STABLE = 1e+8;
	}
	//
	Eigen::VectorXd x_next = x + t*delta_x;
	//
	Eigen::VectorXd feasible = Eigen::VectorXd::Zero(NLINES*NSEGMENTS);
	for(int l=0; l < NLINES; ++l){
		feasible.segment(l*NSEGMENTS, NSEGMENTS) = x_next.segment(l*CL_PROB_SIZE + 1, NSEGMENTS);
	}
	while(!is_less_than_zero(feasible)){
		Rcpp::Rcout << "\rNot feasible, trying to make feasible, t = " << t;
		t *= beta;
		x_next = x + t*delta_x;
		for(int l=0; l < NLINES; ++l){
			feasible.segment(l*NSEGMENTS, NSEGMENTS) = x_next.segment(l*CL_PROB_SIZE + 1, NSEGMENTS);
		}
	}
	//
	double f_x_next = barrier_t*(0.5*(1.0 / NOBS)*x_next.dot(AtA*x_next) - (1.0 / NOBS)*Atb.dot(x_next));
	for(int l=0; l < NLINES; ++l){
		f_x_next += barrier_t*0.5*LAMBDA_G*(1.0 / (NLINES*NGENES))*x_next.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES).squaredNorm() + barrier_t*0.5*LAMBDA_Smooth*(1.0 / (NLINES*NSEGMENTS))*(x_next.segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).dot(SmoothMat*(x_next.segment(l*CL_PROB_SIZE + 1, NSEGMENTS))) - ((-x_next).segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).array().log().matrix().sum();
	}
	f_x_next += barrier_t*0.5*LAMBDA_STABLE*(1.0 / (NGENES))*x.tail(NGENES).squaredNorm();
	while(f_x_next > f_x + alpha*t*(grad_f.dot(delta_x))){
		//
		Rcpp::Rcout << "Backtracking..." << "t = " << t << std::endl;
		//
		t *= beta;
		x_next = x + t*delta_x;
		//
		f_x_next = barrier_t*(0.5*(1.0 / NOBS)*x_next.dot(AtA*x_next) - (1.0 / NOBS)*Atb.dot(x_next));
		for(int l=0; l < NLINES; ++l){
			f_x_next += barrier_t*0.5*LAMBDA_G*(1.0 / (NLINES*NGENES))*x_next.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES).squaredNorm() + barrier_t*0.5*LAMBDA_Smooth*(1.0 / (NLINES*NSEGMENTS))*(x_next.segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).dot(SmoothMat*(x_next.segment(l*CL_PROB_SIZE + 1, NSEGMENTS))) - ((-x_next).segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).array().log().matrix().sum();
		}
		f_x_next += barrier_t*0.5*LAMBDA_STABLE*(1.0 / (NGENES))*x.tail(NGENES).squaredNorm();
		//
		n_steps += 1;
		//
	}
	//
	return(t);
}

Eigen::VectorXd Newton(Eigen::SparseMatrix<double> & AtA, Eigen::VectorXd & Atb, Eigen::SparseMatrix<double> & H, Eigen::VectorXd & g, Eigen::VectorXd & x, Eigen::MatrixXd & SmoothMat,
						double tol, double alpha, double beta, double f_x, int NSEGMENTS, int NLINES, int NGENES, int NOBS, int CL_PROB_SIZE, double barrier_t, double LAMBDA_G, double LAMBDA_Smooth){
	//
	double LAMBDA_STABLE = 1e-8;
	if(NLINES == 1){
		LAMBDA_STABLE = 1e+8;
	}
	double nt_dec = 1;
	double t = 1;
	int n_steps = 1;
	int solver_max_iter = 10;
	Eigen::VectorXd x_nt = Eigen::VectorXd::Zero(x.size());
	Eigen::VectorXd guess = Eigen::VectorXd::Zero(x.size());
	//
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	solver.analyzePattern(H);
	if(solver.info() == Eigen::Success){
		Rcpp::Rcout << "Pattern analysis successful" << std::endl;
	}
	//
	while(TRUE){
		//
		solver.factorize(H);
		if(solver.info() != Eigen::Success){
			//
			Rcpp::Rcout << "\nFactorization unsuccessful" << std::endl;
			//
		}
		if(n_steps > 1 && x_nt.dot(g) <= 0 && n_steps < 20){
			solver_max_iter = 10;
			guess = -1*(x_nt.dot(g))/(x_nt.dot(H*x_nt))*x_nt;
		} else {
			guess = Eigen::VectorXd::Zero(x.size());
			solver_max_iter = 50;
		}
		//
		solver.setMaxIterations(solver_max_iter);
		solver.setTolerance(1e-8);
		Eigen::setNbThreads(8);
		x_nt = solver.solveWithGuess(-g, guess);
		//
		nt_dec = -g.dot(x_nt);
		//
		if(nt_dec < 0){
			Rcpp::Rcout << "\n\nBad Newton dec" << std::endl;
			Rcpp::Rcout << nt_dec << std::endl;
		}
		//
		Rcpp::Rcout << "Finished Newton step " << n_steps << " with decrement " << nt_dec << ", t " << t << ", |g| " << g.squaredNorm() <<  ", and e " << solver.error() << "            \r";
		//
		if(g.squaredNorm() < 1e-4 || nt_dec/2 <= tol || t < 1e-9 || (nt_dec/2 <= 1e-2 && (x_nt.squaredNorm() / x_nt.size()) < 1e-6) || (nt_dec/2 <= 1e-2 && n_steps > 100) || ((x_nt.squaredNorm() / x_nt.size()) < 1e-6 && n_steps > 50)){
			//
			Rcpp::Rcout << "Newton's method converged in " << n_steps << " steps...                                                            " << std::endl;
			return(x);
		//
		} else {
			//
			t = backtrack(AtA, Atb, x, x_nt, g, SmoothMat, alpha, beta, f_x, NSEGMENTS, NLINES, NGENES, NOBS, CL_PROB_SIZE, barrier_t, LAMBDA_G, LAMBDA_Smooth);
			//
			// Compute new Hessian
			//
			for(int w=0; w < NLINES*CL_PROB_SIZE; ++w){
				for(Eigen::SparseMatrix<double>::InnerIterator it(H, w); it; ++it){
					//
					if((it.row() == it.col()) && (((it.row() % CL_PROB_SIZE) > 0) && ((it.row() % CL_PROB_SIZE) < NSEGMENTS + 1))){
						it.valueRef() += std::pow((x(it.row()) + t*x_nt(it.row())), -2) - std::pow(x(it.row()), -2);
					}
					//
				}
			}
			//
			H.makeCompressed();
			//
			// Update x
			//
			x += t*x_nt;
			//
			Eigen::VectorXd feasible = Eigen::VectorXd::Zero(NLINES*NSEGMENTS);
			for(int l=0; l < NLINES; ++l){
				feasible.segment(l*NSEGMENTS, NSEGMENTS) = x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS);
			}
			if(!is_less_than_zero(feasible)){
				Rcpp::Rcout << "Starting with infeasible X, bad, bad, bad" << std::endl;
				for(int i=0; i < feasible.size(); ++i){
					if(feasible(i) >= 0){
						Rcpp::Rcout << "infeasible element is " << feasible(i) << " at position " << i << std::endl;
					}
				}
			}
			//
			// Compute new gradient
			//
			Eigen::VectorXd dphi = Eigen::VectorXd::Zero(x.size());
			for(int l=0; l < NLINES; ++l){
				dphi.segment(l*CL_PROB_SIZE + 1, NSEGMENTS) = barrier_t*LAMBDA_Smooth*(1.0 / (NLINES*NSEGMENTS))*SmoothMat*(x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS)) + (-x).segment(l*CL_PROB_SIZE + 1, NSEGMENTS).array().pow(-1).matrix();
				if(LAMBDA_G != 0){
					dphi.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES) = barrier_t*LAMBDA_G*(1.0 / (NLINES*NGENES))*x.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES);
				}
			}
			if(LAMBDA_STABLE != 0){
				dphi.tail(NGENES) = t*LAMBDA_STABLE*(1.0 / (NGENES))*x.tail(NGENES);
			}
			g = barrier_t*(1.0 / NOBS)*(AtA*x - Atb) + dphi;
			//
			// Compute new f_x
			//
			f_x = barrier_t*(0.5*(1.0 / NOBS)*x.dot(AtA*x) - (1.0 / NOBS)*Atb.dot(x));
			for(int l=0; l < NLINES; ++l){
				f_x += barrier_t*0.5*LAMBDA_G*(1.0 / (NLINES*NGENES))*x.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES).squaredNorm() + barrier_t*0.5*LAMBDA_Smooth*(1.0 / (NLINES*NSEGMENTS))*(x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).dot(SmoothMat*(x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS))) - ((-x).segment(1 + l*CL_PROB_SIZE, NSEGMENTS)).array().log().matrix().sum();
			}
			f_x += barrier_t*0.5*LAMBDA_STABLE*(1.0 / (NGENES))*x.tail(NGENES).squaredNorm();
			//
			n_steps += 1;
		}
	}
}

Eigen::VectorXd Barrier(Eigen::SparseMatrix<double> & AtA, Eigen::VectorXd & Atb, Eigen::VectorXd & x, Eigen::MatrixXd & SmoothMat,
							double t, double mu, double tol, double alpha, double beta, int NSEGMENTS, int NLINES, int NGENES, int NOBS, int CL_PROB_SIZE, double LAMBDA_G, double LAMBDA_Smooth){
	
	int n_steps= 1;
	double LAMBDA_STABLE = 1e-8;
	if(NLINES == 1){
		LAMBDA_STABLE = 1e+8;

	}

	while(TRUE){
		//
		// Compute Hessian at x
		//
		Eigen::SparseMatrix<double> H = t*(1.0 / NOBS)*AtA;
		for(int w=0; w < NLINES*CL_PROB_SIZE; ++w){
			for(Eigen::SparseMatrix<double>::InnerIterator it(H, w); it; ++it){
				//
				if((it.row() == it.col()) && (((it.row() % CL_PROB_SIZE) > 0) && ((it.row() % CL_PROB_SIZE) < NSEGMENTS + 1))){
					//
					it.valueRef() += std::pow(x(it.row()), -2);
					//
				}
				//
				if(LAMBDA_G != 0){
					//
					// Regularize G solutions
					//
					if(it.row() == it.col() && ((it.row() % CL_PROB_SIZE) > NSEGMENTS)){
						it.valueRef() += t*LAMBDA_G*(1.0 / (NLINES*NGENES));
					}
				}
			}
		}
		for(int w=NLINES*CL_PROB_SIZE; w < H.outerSize(); ++w){
			for(Eigen::SparseMatrix<double>::InnerIterator it(H, w); it; ++it){
				//
				// Make gene average terms numerically stable
				//
				if(it.row() == it.col()){
					it.valueRef() += t*LAMBDA_STABLE*(1.0 / (NGENES));
				}
			}
		}
		//
		// Add smoothness penalty
		//
		if(LAMBDA_Smooth != 0){
			for(int w=0; w < SmoothMat.rows(); ++w){
				for(int z=0; z < SmoothMat.cols(); ++z){
					for(int l=0; l < NLINES; ++l){
						H.coeffRef(w + l*CL_PROB_SIZE + 1, z + l*CL_PROB_SIZE + 1) += t*(1.0 / (NLINES*NSEGMENTS))*LAMBDA_Smooth*SmoothMat(w, z);
					}
				}
			}
		}
		//
		H.makeCompressed();
		//
		// Compute gradient vector at x
		//
		Eigen::VectorXd dphi = Eigen::VectorXd::Zero(x.size());
		for(int l=0; l < NLINES; ++l){
			dphi.segment(l*CL_PROB_SIZE + 1, NSEGMENTS) = t*LAMBDA_Smooth*(1.0 / (NLINES*NSEGMENTS))*SmoothMat*(x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS)) + (-x).segment(l*CL_PROB_SIZE + 1, NSEGMENTS).array().pow(-1).matrix();
			if(LAMBDA_G != 0){
				dphi.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES) = t*LAMBDA_G*(1.0 / (NLINES*NGENES))*x.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES);
			}
		}
		if(LAMBDA_STABLE != 0){
			dphi.tail(NGENES) = t*LAMBDA_STABLE*(1.0 / (NGENES))*x.tail(NGENES);
		}
		Eigen::VectorXd g = t*(1.0 / NOBS)*(AtA*x - Atb) + dphi;
		//
		double f_x = t*(0.5*(1.0 / NOBS)*x.dot(AtA*x) - (1.0 / NOBS)*Atb.dot(x));
		for(int l=0; l < NLINES; ++l){
			f_x += t*0.5*LAMBDA_G*(1.0 / (NLINES*NGENES))*x.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES).squaredNorm() + t*0.5*LAMBDA_Smooth*(1.0 / (NLINES*NSEGMENTS))*(x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).dot(SmoothMat*(x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS))) - ((-x).segment(l*CL_PROB_SIZE + 1, NSEGMENTS)).array().log().matrix().sum();
		}
		f_x += t*0.5*LAMBDA_STABLE*(1.0 / (NGENES))*x.tail(NGENES).squaredNorm();
		//
		Rcpp::Rcout << "\nStarting Newton loop..." << std::endl;
		x = Newton(AtA, Atb, H, g, x, SmoothMat, 1e-2, alpha, beta, f_x, NSEGMENTS, NLINES, NGENES, NOBS, CL_PROB_SIZE, t, LAMBDA_G, LAMBDA_Smooth);
		//
		if(double(NLINES*NSEGMENTS)/t < tol){
			Rcpp::Rcout << "Barrier method converged in " << n_steps << " steps...                   \n";
			return(x);
		} else {
			t *= mu;
			n_steps += 1;
		}
	}
}

double backtrack2(Eigen::SparseMatrix<double> & AtA, Eigen::VectorXd & Atb, Eigen::VectorXd & x, Eigen::VectorXd & delta_x, Eigen::VectorXd & grad_f, 
					double alpha, double beta, double f_x, int NGUIDES, int NOBS, double barrier_t, double LAMBDA_Off){
	//
	double t = 1;
	int n_steps = 1;
	//
	Eigen::VectorXd x_next = x + t*delta_x;
	//
	Eigen::VectorXd feasible = x_next.segment(NGUIDES, NGUIDES);
	while(!is_between_one_and_zero(feasible)){
		t *= beta;
		x_next = x + t*delta_x;
		feasible = x_next.segment(NGUIDES, NGUIDES);
	}
	//
	double f_x_next = barrier_t*(0.5*(1.0 / NOBS)*x_next.dot(AtA*x_next) - (1.0 / NOBS)*Atb.dot(x_next) + 0.5*LAMBDA_Off*(1.0 / NGUIDES)*(x_next.segment(0, NGUIDES).squaredNorm())) - (x_next.segment(NGUIDES, NGUIDES)).array().log().matrix().sum() - (1 - x_next.segment(NGUIDES, NGUIDES).array()).log().matrix().sum();
	while(f_x_next > f_x + alpha*t*(grad_f.dot(delta_x))){
		t *= beta;
		x_next = x + t*delta_x;
		f_x_next = barrier_t*(0.5*(1.0 / NOBS)*x_next.dot(AtA*x_next) - (1.0 / NOBS)*Atb.dot(x_next) + 0.5*LAMBDA_Off*(1.0 / NGUIDES)*(x_next.segment(0, NGUIDES).squaredNorm())) - (x_next.segment(NGUIDES, NGUIDES)).array().log().matrix().sum() - (1 - x_next.segment(NGUIDES, NGUIDES).array()).log().matrix().sum();
		n_steps += 1;
	}
	//
	return(t);
}

Eigen::VectorXd Newton2(Eigen::SparseMatrix<double> & AtA, Eigen::VectorXd & Atb, Eigen::SparseMatrix<double> & H, Eigen::VectorXd & g, Eigen::VectorXd & x,
						double tol, double alpha, double beta, double f_x, int NGUIDES, int NOBS, double barrier_t, double LAMBDA_Off){
	//
	double nt_dec = 1;
	double t = 1;
	int n_steps = 1;
	//
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver;
	solver.analyzePattern(H);
	//
	while(TRUE){
		solver.factorize(H);
		Eigen::VectorXd x_nt = solver.solve(-g);
		nt_dec = g.dot(-x_nt);
		//
		Rcpp::Rcout << "Finished Newton step " << n_steps << " with decrement " << nt_dec << " and t " << t <<  "      \r";
		if(nt_dec/2 <= tol || t < 1e-9){
			Rcpp::Rcout << "Newton's method converged in " << n_steps << " steps...                         " << std::endl;
			return(x);
		} else {
			t = backtrack2(AtA, Atb, x, x_nt, g, alpha, beta, f_x, NGUIDES, NOBS, barrier_t, LAMBDA_Off);
			//
			// Compute new Hessian
			//
			Eigen::VectorXd invx2 = x.segment(NGUIDES, NGUIDES).array().pow(-2).matrix() + (1 - x.segment(NGUIDES, NGUIDES).array()).pow(-2).matrix();
			Eigen::VectorXd invx2_next = (x + t*x_nt).segment(NGUIDES, NGUIDES).array().pow(-2).matrix() + (1 - (x + t*x_nt).segment(NGUIDES, NGUIDES).array()).pow(-2).matrix();
			for(int w=NGUIDES; w < 2*NGUIDES; ++w){
				for(Eigen::SparseMatrix<double>::InnerIterator it(H, w); it; ++it){
					if(it.row() == it.col()){
						it.valueRef() += invx2_next(it.row() - NGUIDES) - invx2(it.row() - NGUIDES);
					}
				}
			}
			//
			// Update x
			//
			x += t*x_nt;
			//
			// Compute new gradient
			//
			Eigen::VectorXd dphi = Eigen::VectorXd::Zero(x.size());
			dphi.segment(NGUIDES, NGUIDES) = (1 - x.segment(NGUIDES, NGUIDES).array()).pow(-1).matrix() - x.segment(NGUIDES, NGUIDES).array().pow(-1).matrix();
			dphi.segment(0, NGUIDES) = barrier_t*LAMBDA_Off*(1.0 / NGUIDES)*x.segment(0, NGUIDES);
			g = barrier_t*(1.0 / NOBS)*(AtA*x - Atb) + dphi;
			//
			// Compute new f_x
			//
			f_x = barrier_t*(0.5*(1.0 / NOBS)*x.dot(AtA*x) - (1.0 / NOBS)*Atb.dot(x) + 0.5*LAMBDA_Off*(1.0 / NGUIDES)*(x.segment(0, NGUIDES).squaredNorm())) - (x.segment(NGUIDES, NGUIDES)).array().log().matrix().sum() - (1 - x.segment(NGUIDES, NGUIDES).array()).log().matrix().sum();
			//
			n_steps += 1;
		}
	}
}

Eigen::VectorXd Barrier2(Eigen::SparseMatrix<double> & AtA, Eigen::VectorXd & Atb, Eigen::VectorXd & x, 
							double t, double mu, double tol, double alpha, double beta, int NGUIDES, int NOBS, double LAMBDA_Off){
	int n_steps= 1;
	while(TRUE){
		//
		// Compute Hessian at x
		//
		Eigen::VectorXd invx2 = x.segment(NGUIDES, NGUIDES).array().pow(-2).matrix() + (1 - x.segment(NGUIDES, NGUIDES).array()).pow(-2).matrix();
		Eigen::SparseMatrix<double> H = t*(1.0 / NOBS)*AtA;
		for(int w=0; w < H.outerSize(); ++w){
			for(Eigen::SparseMatrix<double>::InnerIterator it(H, w); it; ++it){
				if(it.row() == it.col()){
					if(it.row() >= NGUIDES){
						it.valueRef() += invx2(it.row() - NGUIDES);
					} else {
						// add regularization on guide offset
						it.valueRef() += t*(1.0 / NGUIDES)*LAMBDA_Off;
					}
				}
			}
		}
		//
		// Compute gradient vector at x
		//
		Eigen::VectorXd dphi = Eigen::VectorXd::Zero(x.size());
		dphi.segment(NGUIDES, NGUIDES) = (1 - x.segment(NGUIDES, NGUIDES).array()).pow(-1).matrix() - x.segment(NGUIDES, NGUIDES).array().pow(-1).matrix();
		dphi.segment(0, NGUIDES) = t*LAMBDA_Off*(1.0 / NGUIDES)*x.segment(0, NGUIDES);
		Eigen::VectorXd g = t*(1.0 / NOBS)*(AtA*x - Atb) + dphi;
		//
		double f_x = t*(0.5*(1.0 / NOBS)*x.dot(AtA*x) - (1.0 / NOBS)*Atb.dot(x) + 0.5*LAMBDA_Off*(1.0 / NGUIDES)*(x.segment(0, NGUIDES).squaredNorm())) - (x.segment(NGUIDES, NGUIDES)).array().log().matrix().sum() - (1 - x.segment(NGUIDES, NGUIDES).array()).log().matrix().sum();
		//
		x = Newton2(AtA, Atb, H, g, x, tol, alpha, beta, f_x, NGUIDES, NOBS, t, LAMBDA_Off);
		//
		if(double(2*NGUIDES)/t < tol){
			Rcpp::Rcout << "Barrier method converged in " << n_steps << " steps...              \n";
			return(x);
		} else {
			t *= mu;
			n_steps += 1;
		}
	}

}

Eigen::SparseMatrix<double> makeValidationSet(Eigen::MatrixXd & D, Eigen::SparseMatrix<double> & M, double VAL_FRAC=0.1, int NLINES=1){
	//
	std::default_random_engine generator;
	std::uniform_int_distribution<int> guide_distribution(0, std::pow(VAL_FRAC, -1) - 1);
	//
	typedef Eigen::Triplet<double> T;
	std::vector<T> valtripletList;
	//
	Eigen::VectorXd guides_in_val;
	//
	Eigen::VectorXd guides_per_gene = M.transpose()*(Eigen::VectorXd::Ones(M.rows()));
	Eigen::VectorXd genes_per_guide = M*(Eigen::VectorXd::Ones(M.cols()));
	//
	Eigen::MatrixXd gene_guides_in_val;
	//
	int is_valid = 0;
	//
	int n_trys = 1;
	//
	Eigen::SparseMatrix<double> Val(D.rows(), D.cols());
	//
	while(!is_valid){
		//
		// Clear validation set
		//
		valtripletList.clear();
		//
		Eigen::VectorXd guides_visited = Eigen::VectorXd::Zero(D.rows());
		//
		for(int i=0; i < M.outerSize(); ++i){
			if(guides_per_gene(i) < 2){
				for(Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
					guides_visited(it.row()) = 1;
				}
			}
		}
		//
		for(int i=0; i < M.outerSize(); ++i){
			//
			int nguides_left = guides_per_gene(i);
			//
			int add_from_this_gene = 1;
			//
			if(nguides_left < 2){
				add_from_this_gene = 0;
				for(Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
					guides_visited(it.row()) = 1;
				}
			}
			//
			if(add_from_this_gene){
				for(int l=0; l < D.cols(); ++l){
					int guide_to_omit = guide_distribution(generator);
					if(guide_to_omit < nguides_left){
						int guide_index = 0;
						for(Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
							if((guide_to_omit == guide_index) && (guides_visited(it.row()) == 0)){
								valtripletList.push_back(T(it.row(), l, 1));
							}
							guide_index += 1;
						}
					}
				}
				for(Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it){
					guides_visited(it.row()) = 1;
				}
			}
		}
		//
		Val = Eigen::SparseMatrix<double>(D.rows(), D.cols());
		Val.setFromTriplets(valtripletList.begin(), valtripletList.end());
		//
		for(int i=0; i < Val.outerSize(); ++i){
			for(Eigen::SparseMatrix<double>::InnerIterator it(Val, i); it; ++it){
				if(it.value() > 1){
					it.valueRef() = 1;
				}
			}
		}
		//
		// Check that there are no guides that are entirely in the validation set
		//
		int is_guide_valid = 1;
		guides_in_val = Val*Eigen::VectorXd::Ones(D.cols());
		for(int j=0; j < guides_in_val.size(); ++j){
			if(guides_in_val(j) == D.cols()){
				is_guide_valid = 0;
			}
		}
		//
		// Check that no gene / cell line pairs are entirely in validation set
		//
		int is_gene_valid = 1;
		gene_guides_in_val = (M.transpose())*Val;
		int num_invalid = 0;
		for(int k=0; k < gene_guides_in_val.rows(); ++k){
			for(int l=0; l < gene_guides_in_val.cols(); ++l){
				if(gene_guides_in_val(k, l) == guides_per_gene(k) && (guides_per_gene(k) != 0)){
					is_gene_valid = 0;
					num_invalid += 1;
				}
			}
		}
		//
		if(NLINES < 5){
			is_guide_valid = 1;
		}
		//
		is_valid = is_guide_valid*is_gene_valid;
		//
		n_trys = n_trys + 1;
		//
	}
	Rcpp::Rcout << "Made validation set in " << n_trys << " attempts" << std::endl;
	Rcpp::Rcout << "Validation set fraction is " << double(Val.nonZeros()) / (Val.rows()*Val.cols()) << " and desired is " << VAL_FRAC << std::endl;
	return(Val);
}


Eigen::VectorXd findRedundantGenes(Eigen::SparseMatrix<double> & MtM){
	//
	typedef Eigen::Triplet<double> T;
	//
	// Check to see that M is full column rank
	//
	MtM.makeCompressed();
	Eigen::SparseMatrix<double> BadGenes(MtM.cols(), MtM.cols());
	std::vector<T> BadtripletList;
	BadtripletList.reserve(MtM.cols());
	//
	Rcpp::Rcout << "Finding redundant columns..." << std::endl;
	int num_bad_genes = 0;
	for(int k=0; k < MtM.outerSize(); ++k){
		int num_in_col = 0;
		for(Eigen::SparseMatrix<double>::InnerIterator it(MtM, k); it; ++it){
			if(it.value() != 0){
				num_in_col += 1;
			}
		}
		Eigen::VectorXd MtMcolumn = Eigen::VectorXd(num_in_col);
		int counter = 0;
		for(Eigen::SparseMatrix<double>::InnerIterator it(MtM, k); it; ++it){
			if(it.value() != 0){
				MtMcolumn(counter) = it.row();
				counter += 1;
			}
		}
		for(int l=0; l < MtM.outerSize(); ++l){
			if(l != k){
				int num_in_col2 = 0;
				for(Eigen::SparseMatrix<double>::InnerIterator it2(MtM, l); it2; ++it2){
					if(it2.value() != 0){
						num_in_col2 += 1;
					}
				}
				if(num_in_col2 == num_in_col){
					int equals_other_column = 1;
					int counter2 = 0;
					for(Eigen::SparseMatrix<double>::InnerIterator it2(MtM, l); it2; ++it2){
						if(it2.value() != 0){
							if(it2.row() != MtMcolumn(counter2)){
								equals_other_column = 0;
							}
							counter2 += 1;
						}
					}
					if(equals_other_column){
						num_bad_genes += 1;
						BadtripletList.push_back(T(k, l, 1));
					}
				}
			}
		}
	}
	//
	BadGenes.setFromTriplets(BadtripletList.begin(), BadtripletList.end());
	//
	for(int l=0; l < BadGenes.outerSize(); ++l){
		for(Eigen::SparseMatrix<double>::InnerIterator it2(BadGenes, l); it2; ++it2){
			if(it2.valueRef() != 0){
				it2.valueRef() = 1;
			}
		}
	}
	//
	Rcpp::Rcout << "There are " << (BadGenes*Eigen::VectorXd::Ones(BadGenes.cols())).sum() / 2 << " redundant genes total " << std::endl;
	//
	// Find the redundant genes for removal
	//
	Eigen::VectorXd genes_to_remove = Eigen::VectorXd::Ones(BadGenes.cols());
	for(int l=0; l < BadGenes.outerSize(); ++l){
		if(genes_to_remove(l) != 0){
			genes_to_remove(l) = -1;
		}
		for(Eigen::SparseMatrix<double>::InnerIterator it2(BadGenes, l); it2; ++it2){
			if((it2.value() != 0) && (genes_to_remove(it2.row()) != -1)){
				genes_to_remove(it2.row()) = 0;
			}
		}
	}
	genes_to_remove = (genes_to_remove.array().abs() - 1).abs().matrix();
	Rcpp::Rcout << "Will remove " << genes_to_remove.sum() << " redundant genes from model..." << std::endl;
	//
	return(genes_to_remove);
}

Eigen::MatrixXd removeRedundantGenes(Eigen::MatrixXd & G, Eigen::VectorXd & genes_to_remove){
	//
	Eigen::MatrixXd Gnew = Eigen::MatrixXd::Zero(G.rows() - genes_to_remove.sum(), G.cols());
	int counter = 0;
	for(int i=0; i < G.rows(); ++i){
		if(!genes_to_remove(i)){
			Gnew.row(counter) = G.row(i);
			counter += 1;
		}
	}
	//
	return(Gnew);
	//
}

Eigen::SparseMatrix<double, Eigen::RowMajor> removeRedundantMappings(Eigen::SparseMatrix<double, Eigen::RowMajor> & M,
																		Eigen::VectorXd & genes_to_remove){
	//
	typedef Eigen::Triplet<double> T;
	//
	Eigen::SparseMatrix<double, Eigen::RowMajor> Mnew(M.rows(), M.cols() - genes_to_remove.sum());
	std::vector<T> MnewtripletList;
	MnewtripletList.reserve(M.nonZeros());
	//
	for(int k=0; k < M.outerSize(); ++k){
		//
		for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(M, k); it; ++it){
			//
			if(!genes_to_remove(it.col())){
				//
				MnewtripletList.push_back(T(it.row(), it.col() - genes_to_remove.head(it.col()+1).sum(), 1));
				//
			}
			//
		}
		//
	}
	//
	Mnew.setFromTriplets(MnewtripletList.begin(), MnewtripletList.end());
	//
	return(Mnew);
	//
}

void writeTextToLog(std::string & text, std::string & name){
	std::ofstream log_file(name, std::ios_base::out | std::ios_base::app);
	log_file << text << std::endl;
}

Eigen::RowVectorXd makeQuantiles(Eigen::VectorXd cn_vec, int NSEGMENTS, int min_scale, int max_scale){
	//
	std::sort(cn_vec.data(), cn_vec.data() + cn_vec.size());
	//
	double min_resolution = double(cn_vec.maxCoeff() - cn_vec.minCoeff()) / double(min_scale);
	double max_resolution = double(cn_vec.maxCoeff() - cn_vec.minCoeff()) / double(max_scale);
	Rcpp::Rcout << "Min resolution is " << min_resolution << " copies..." << std::endl;
	Rcpp::Rcout << "Max resolution is " << max_resolution << " copies..." << std::endl;
	Rcpp::Rcout << "Max entry " << cn_vec.maxCoeff() << " last entry " << cn_vec(cn_vec.size()-1) << std::endl;
	Rcpp::Rcout << "Min entry " << cn_vec.minCoeff() << " first entry " << cn_vec(0) << std::endl;
	//
	int n_quantiles = NSEGMENTS;
	Eigen::RowVectorXd quantile_vec = Eigen::RowVectorXd::Zero(2*NSEGMENTS);
	for(int i=0; i < NSEGMENTS; ++i){
		quantile_vec(2*i) = cn_vec(std::floor(double(i*cn_vec.size())/n_quantiles));
		quantile_vec(2*i + 1) = cn_vec(std::floor(double((i+1)*cn_vec.size())/n_quantiles));
	}
	quantile_vec(0) = 0;
	quantile_vec(quantile_vec.size() - 1) = cn_vec.maxCoeff() + 1;
	//
	double threshold = 1e+6;
	while((threshold > min_resolution) && (n_quantiles < cn_vec.size())){
		//
		double max_quantile_span = 0;
		for(int i=0; i < n_quantiles; ++i){
			int right_index = std::floor(double((i+1)*cn_vec.size())/n_quantiles);
			int left_index = std::floor(double(i*cn_vec.size())/n_quantiles);
			int vec_last_entry = cn_vec.size() - 1;
			double quantile_span = cn_vec(std::min(right_index, vec_last_entry)) - cn_vec(std::min(left_index, vec_last_entry));
			if(quantile_span > max_quantile_span){
				max_quantile_span = quantile_span;
			}
		}
		//
		threshold = max_quantile_span;
		//
		n_quantiles *= 1.2;
	//
	}
	Rcpp::Rcout << n_quantiles << " quantiles with max span " << threshold << std::endl;
	//
	Eigen::VectorXd breaks = Eigen::VectorXd::Zero(n_quantiles + 1);
	for(int i=0; i < breaks.size(); ++i){
		int left_index = std::floor(double(i*cn_vec.size())/n_quantiles);
		breaks(i) = cn_vec(left_index);
	}
	breaks(0) = 0;
	breaks(breaks.size() - 1) = cn_vec.maxCoeff();
	//
	Eigen::VectorXd spans = breaks.tail(n_quantiles) - breaks.head(n_quantiles);
	//
	double distance_threshold = spans.maxCoeff();
	//
	int num_clusters = 1;
	while(TRUE){
		//
		num_clusters = 1;
		for(int i=0; i < spans.size(); ++i){
			if((spans(i) > distance_threshold) && (breaks(i) > 0)){
				num_clusters += 1;
			}
		}
		if(num_clusters == NSEGMENTS){
			break;
		} else if(num_clusters < NSEGMENTS){
			distance_threshold *= 0.95;
		} else {
			distance_threshold *= 1.01;
		}
		Rcpp::Rcout << "\r" << distance_threshold;
	//
	}
	Rcpp::Rcout << std::endl;
	Rcpp::Rcout << num_clusters << " clusters makes " << NSEGMENTS << " segments with distance " << distance_threshold << "..." << std::endl;
	//
	Rcpp::Rcout << "spans size is " << spans.size() << " and breaks size is " << breaks.size() << std::endl;
	int counter = 0;
	for(int i=0; i < spans.size(); ++i){
		if((spans(i) > distance_threshold) && (breaks(i) > 0)){
			//
			if(counter < NSEGMENTS - 1){
				quantile_vec(2*counter + 1) = breaks(i);
				quantile_vec(2*(counter + 1)) = breaks(i);
			} else {
				quantile_vec(2*counter + 1) = breaks(i);
			}
			//
			counter += 1;
		}
		//
	}
	//
	Rcpp::Rcout << quantile_vec << std::endl;
	return(quantile_vec);
//
}


// Pass matrices by reference so that they are modified
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List fit_ceres(NumericMatrix & rD,
			  NumericMatrix & rQ, 
			  NumericMatrix & rM,
			  NumericVector & rColCl, 
			  NumericMatrix & rG, 
			  NumericMatrix & rC, NumericVector & rTox,
			  NumericMatrix & rQuantileMat,
			  double LAMBDA_G, double LAMBDA_Off, 
			  double LAMBDA_Smooth, int NSEGMENTS,
			  int MAKE_VALIDATION_SET,
			  String log_file_suffix = "",
			  String log_file_dir = "log",
			  int fit_efficacy = 1)
{

	// Set random seed
	std::default_random_engine generator;
	generator.seed(time(NULL));
	//
	// Define number of genes, loci, cell lines
	int NGENES = rG.nrow();
	const int NGUIDES = rD.nrow();
	const int NLOCI = rC.nrow();
	const int NLINES = rG.ncol();
	const int NREPS = rD.ncol();
	const double VAL_FRAC = 0.1;
	const double ALPHA = 0.15;
	const double BETA = 0.5;
	const double TOL = 1e-4;
	const double MU = 10;
	const double t0 = ((NLINES*NSEGMENTS) / MU);
	const double t02 = (2*NGUIDES / MU);
	//
	// Define Triplet class
	typedef Eigen::Triplet<double> T;
	std::vector<T> MtripletList;
	MtripletList.reserve(rM.nrow());
	std::vector<T> QindtripletList;
	QindtripletList.reserve(rQ.nrow());
	//
	// Initialize matrices
	Eigen::VectorXd Q = 0.9999*Eigen::VectorXd::Ones(NGUIDES);
	Eigen::MatrixXd G = Eigen::MatrixXd::Random(NGENES, NLINES);
	Eigen::MatrixXd Tox = -1*(Eigen::MatrixXd::Random(NLINES, NSEGMENTS)).cwiseAbs();
	Eigen::RowVectorXd I = 0.1*Eigen::RowVectorXd::Random(NLINES).cwiseAbs();
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(NLOCI, NLINES);
	Eigen::MatrixXd D = Eigen::MatrixXd::Zero(NGUIDES, NREPS);
	Eigen::VectorXd Off = Eigen::VectorXd::Zero(NGUIDES);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(NGENES);
	Eigen::MatrixXd Dhat;
	Eigen::VectorXd ColCl = Eigen::VectorXd::Zero(NREPS);
	double train_error = 0;
	double val_error = 0;
	double temp_error = 0;
	double loss = 0;
	//
	//
	// M
	Rcpp::Rcout << "Converting M" << std::endl;
	for(int i=0; i < rM.nrow(); ++i){
		MtripletList.push_back(T(rM(i,0), rM(i,1), rM(i,2)));
	}
	Eigen::SparseMatrix<double, Eigen::RowMajor> M(NLOCI, NGENES);
	Eigen::SparseMatrix<double, Eigen::ColMajor> Mcol(NLOCI, NGENES);
	M.setFromTriplets(MtripletList.begin(), MtripletList.end());
	M.makeCompressed();
	Mcol.setFromTriplets(MtripletList.begin(), MtripletList.end());
	//
	// Q
	Rcpp::Rcout << "Converting Q" << std::endl;
	for(int i=0; i < rQ.nrow(); ++i){
		QindtripletList.push_back(T(rQ(i,0), rQ(i,1), 1));
	}
	Eigen::SparseMatrix<double, Eigen::RowMajor> Qind(NGUIDES, NLOCI);
	Eigen::SparseMatrix<double, Eigen::ColMajor> Qindcol(NGUIDES, NLOCI);
	Qind.setFromTriplets(QindtripletList.begin(), QindtripletList.end());
	Qind.makeCompressed();
	Qindcol.setFromTriplets(QindtripletList.begin(), QindtripletList.end());
	//
	// C
	Rcpp::Rcout << "Converting C" << std::endl;
	for(int i=0; i < C.rows(); ++i){
		for(int j=0; j < C.cols(); ++j){
			C(i, j) = rC(i, j);
		}
	}
	//
	// D
	Rcpp::Rcout << "Converting D" << std::endl;
	for(int i=0; i < D.rows(); ++i){
		for(int j=0; j < D.cols(); ++j){
			if(std::isnan(rD(i,j))){
				D(i, j) = 0;
			} else {
				D(i, j) = rD(i, j);
			}
		}
	}
	//
	// ColCl
	for(int i=0; i < ColCl.size(); ++i){
		ColCl(i) = rColCl[i];
	}
	//
	// Perform last-minute initializations
	//
	M = Qind*M;
	Mcol = Qindcol*Mcol;
	Eigen::MatrixXd Ceff = Qind*C;
	for(int k=0; k < M.outerSize(); ++k){
		for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(M, k); it; ++it){
			if(it.value() > 1){
				it.valueRef() = 1;
			}
		}
	}
	for(int k=0; k < Mcol.outerSize(); ++k){
		for(Eigen::SparseMatrix<double>::InnerIterator it(Mcol, k); it; ++it){
			if(it.value() > 1){
				it.valueRef() = 1;
			}
		}
	}
	M.makeCompressed();
	Mcol.makeCompressed();
	//
	Eigen::MatrixXd Mtriplets = Eigen::MatrixXd::Zero(M.nonZeros(), 3);
	int index = 0;
	for(int k=0; k < M.outerSize(); ++k){
		for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(M, k); it; ++it){
			Mtriplets(index, 0) = it.row();
			Mtriplets(index, 1) = it.col();
			Mtriplets(index, 2) = it.value();
		}
	}
	Eigen::MatrixXd CeffTox = Eigen::MatrixXd::Zero(NGUIDES, NLINES);
	Eigen::MatrixXd QuantileMat = Eigen::MatrixXd::Zero(NLINES, 2*NSEGMENTS);
	for(int i=0; i < rQuantileMat.nrow(); ++i){
		for(int j=0; j < rQuantileMat.ncol(); ++j){
			QuantileMat(i, j) = rQuantileMat(i, j);
		}
	}
	//
	Eigen::MatrixXd SmoothMat;
	Eigen::MatrixXd Smoothness;
	if(NSEGMENTS > 1){
		Smoothness = Eigen::MatrixXd::Zero(NSEGMENTS, NSEGMENTS);
		for(int s=0; s < NSEGMENTS-1; ++s){
			Smoothness(s, s) = -1.0;
			Smoothness(s, s+1) = 1.0;
		}
		Smoothness(NSEGMENTS-1, NSEGMENTS-1) = 1;
		SmoothMat = (Smoothness.transpose()*Smoothness);
	} else {
		SmoothMat = Eigen::MatrixXd::Ones(1, 1);
	}
	//
	Eigen::SparseMatrix<double> Val = Eigen::SparseMatrix<double>(NGUIDES, NREPS);

	if(MAKE_VALIDATION_SET){
		Val = makeValidationSet(D, Mcol, VAL_FRAC, NLINES);
	} else {
		std::vector<T> valtripletList;
		valtripletList.push_back(T(0, 0, 0));
		Val.setFromTriplets(valtripletList.begin(), valtripletList.end());
	}
	//
	Eigen::VectorXd reps_per_cl = Eigen::VectorXd::Zero(NLINES);
	for(int l=0; l < NLINES; ++l){
		int num_reps = 0;
		for(int r=0; r < ColCl.size(); ++r){
			if(ColCl(r) == l){
				num_reps += 1;
			}
		}
		reps_per_cl(l) = num_reps;
	}
	//
	double iter = 0;
	//
	// Find redundant genes in M
	//
	Eigen::VectorXd genes_to_remove = findRedundantGenes(Mcol);
	//
	// Remove redundant genes in G and M
	//
	Rcpp::Rcout << "G has " << G.rows() << " rows " << std::endl;
	G = removeRedundantGenes(G, genes_to_remove);
	NGENES = G.rows();
	Rcpp::Rcout << "G has " << G.rows() << " rows " << std::endl; 
	Rcpp::Rcout << "Removed redundant genes..." << std::endl;
	Rcpp::Rcout << "M has " << M.cols() << " cols " << std::endl; 
	Eigen::SparseMatrix<double, Eigen::RowMajor> Mnew = removeRedundantMappings(M, genes_to_remove);
	Rcpp::Rcout << "Removed redundant mappings..." << std::endl;
	Rcpp::Rcout << "M has " << Mnew.cols() << " cols " << std::endl;
	//
	// Initialize log file
	//
	std::string log_file_suffix_cpp = log_file_suffix;
	std::string log_file_dir_cpp = log_file_dir;
	if(log_file_suffix_cpp != ""){
		log_file_suffix_cpp = "_" + log_file_suffix_cpp;
	}
	std::string logfileName = log_file_dir_cpp + "/logfile" + log_file_suffix_cpp + ".txt";
	Rcpp::Rcout << logfileName << std::endl;
	std::string colnames = std::string("loss\tloss_g\tloss_s\ttrain_err\ttest_err");
	//
	// Write colnames to logfile
	//
	writeTextToLog(colnames, logfileName);
	//
	int MAXITER;
	if(NLINES < 5){
		MAXITER = 1;
	} else {
		MAXITER = 5;
	}
	//
	while(iter < MAXITER){
		//
		Eigen::MatrixXd Val_Set = Val.cwiseProduct(Eigen::MatrixXd::Ones(NGUIDES, NREPS));
		//
		// Report error
		//
		Dhat = Off.rowwise().replicate(Ceff.cols()) + I.colwise().replicate(Q.size()) + Q.asDiagonal()*(Mnew*G  + CeffTox);
		train_error = 0;
		val_error = 0;
		temp_error = 0;
		//
		for(int c=0; c < NREPS; ++c){
			temp_error = ((1.0/std::sqrt(reps_per_cl(ColCl(c))))*(D.col(c) - Dhat.col(ColCl(c)))).cwiseProduct(Val.col(c)).squaredNorm();
			val_error +=  temp_error / Val_Set.col(c).sum();
			train_error += (((1.0/std::sqrt(reps_per_cl(ColCl(c))))*(D.col(c) - Dhat.col(ColCl(c)))).squaredNorm() - temp_error);
		}
		val_error = val_error / NLINES;
		train_error = train_error / (NGUIDES*NLINES);
		Rcpp::Rcout << "\nTraining error is " << train_error << std::endl;
		Rcpp::Rcout << "Test error is " << val_error << std::endl;
		Rcpp::Rcout << "Smoothness is " << (Smoothness*(Tox.transpose())).squaredNorm() << std::endl;
		Rcpp::Rcout << "Gene sol size is " << G.squaredNorm() << std::endl;
		//
		// Log error
		//
		loss = train_error + 0.5*LAMBDA_Smooth*(1.0 / NSEGMENTS)*(Smoothness*(Tox.transpose())).squaredNorm() + 0.5*LAMBDA_G*(1.0 / (G.rows()*G.cols()))*G.squaredNorm();
		std::string lossMessage = std::to_string(loss);
		std::string lossGMessage = "\t" + std::to_string(0.5*LAMBDA_G*(1.0 / (G.rows()*G.cols()))*G.squaredNorm());
		std::string lossSMessage = "\t" + std::to_string(0.5*LAMBDA_Smooth*(1.0 / NSEGMENTS)*(Smoothness*(Tox.transpose())).squaredNorm());
		std::string trainingErrorMessage = "\t" + std::to_string(train_error);
		std::string validationErrorMessage = "\t" + std::to_string(val_error);
		//
		std::string message = lossMessage + lossGMessage + lossSMessage + trainingErrorMessage + validationErrorMessage;
		writeTextToLog(message, logfileName);
		//
		// Solve sparse least squares problem for gene solutions and cutting toxicity
		//
		Rcpp::Rcout << "\nSolving for gene solutions and toxicities... " << std::endl;
		//
		// Update M
		//
		Eigen::SparseMatrix<double, Eigen::RowMajor> Mtemp = Q.asDiagonal()*Mnew;
		//
		// Update Ceff
		//
		Eigen::VectorXd Ceffvec(Eigen::Map<Eigen::VectorXd>(Ceff.data(), NGUIDES*NLINES));
		//
		int VAL_SIZE = Val_Set.colwise().sum().sum();
		int val_sum = 0;
		//
		// Instantiate big sparse matrix
		//
		int CL_PROB_SIZE = NGENES + NSEGMENTS + 1;
		Eigen::SparseMatrix<double> A(NGUIDES*NREPS - VAL_SIZE, NLINES*(CL_PROB_SIZE) + NGENES);
		std::vector<T> AtripletList;
		AtripletList.reserve(2*NREPS*Mtriplets.size() + NREPS*NSEGMENTS + NREPS);
		//
		for(int index=0; index < NGUIDES*NREPS; ++index){
			// Set guide / cell line counter
			int gl = (index % NGUIDES) + NGUIDES*ColCl(index / NGUIDES);
			//
			if(Val_Set(index % NGUIDES, index / NGUIDES) == 0){
				// add intercept term
				AtripletList.push_back(T(val_sum, (gl / NGUIDES)*(CL_PROB_SIZE), 1));
				// add effective cuts term for cell lines
				for(int s=0; s < NSEGMENTS; ++s){
					if((Ceffvec(gl) >= QuantileMat(gl / NGUIDES, 2*s)) && (Ceffvec(gl) < QuantileMat(gl / NGUIDES, 2*s + 1)) && (Ceffvec(gl) != 0)){
						for(int c=0; c < s; ++c){
							AtripletList.push_back(T(val_sum, (gl / NGUIDES)*(CL_PROB_SIZE) + c + 1, Q(gl % NGUIDES)*(QuantileMat(gl / NGUIDES, 2*c + 1) - QuantileMat(gl / NGUIDES, 2*c))));
						}
						AtripletList.push_back(T(val_sum, (gl / NGUIDES)*(CL_PROB_SIZE) + s + 1, Q(gl % NGUIDES)*(Ceffvec(gl) - QuantileMat(gl / NGUIDES, 2*s))));
					}
				}
				// add guide-to-gene mapping terms
				for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Mtemp, gl % NGUIDES); it; ++it){
					if(it.value() != 0){
						// One for variable gene component
						AtripletList.push_back(T(val_sum, it.col() + (gl / NGUIDES)*(CL_PROB_SIZE) + NSEGMENTS + 1, it.value()));
						// One for gene-wide mean across cell lines
						AtripletList.push_back(T(val_sum, it.col() + NLINES*(CL_PROB_SIZE), it.value()));
					}
				}
				// increment val_sum
				val_sum += 1;
			}
		}
		//
		A.setFromTriplets(AtripletList.begin(), AtripletList.end());
		A.makeCompressed();
		Rcpp::Rcout << "Instantiated model matrix..." << std::endl;
		//
		// Weight replicates so that cell lines have equal weight in minimizing error
		//
		Eigen::VectorXd rep_weight = Eigen::VectorXd::Zero(A.rows());
		val_sum = 0;
		for(int gl=0; gl < NGUIDES*NREPS; ++gl){
			if(Val_Set(gl % NGUIDES, gl / NGUIDES) == 0){
				//
				rep_weight(val_sum) = std::sqrt(1.0 / reps_per_cl(ColCl(gl / NGUIDES)));
				// Increment val_sum
				val_sum += 1;
			}
		}
		//
		// Reweight A
		//
		A = rep_weight.asDiagonal()*A;
		//
		// Compute A'A
		//
		Eigen::SparseMatrix<double> AtA = A.transpose()*A;
		//
		// Set target vector
		//
		Eigen::VectorXd b(Eigen::Map<Eigen::VectorXd>(D.data(), D.rows()*D.cols()));
		Eigen::VectorXd new_b = Eigen::VectorXd::Zero(NGUIDES*NREPS - VAL_SIZE);
		b -= Off.colwise().replicate(NREPS);
		val_sum = 0;
		for(int gl=0; gl < NGUIDES*NREPS; ++gl){
			if(Val_Set(gl % NGUIDES, gl / NGUIDES) == 0){
				//
				new_b(val_sum) = b(gl);
				// Increment val_sum
				val_sum += 1;
			}
		}
		//
		// Reweight new_b
		//
		new_b = rep_weight.asDiagonal()*new_b;
		//
		// Set target vector
		//
		Eigen::VectorXd Atb = A.transpose()*new_b;
		//
		// Instantiate gene solution vector, requiring that it is feasible, i.e. all toxicity params are < 0
		//
		Eigen::VectorXd x = Eigen::VectorXd::Zero(NLINES*(CL_PROB_SIZE) + NGENES);
		for(int l=0; l < NLINES; ++l){
			x(l*CL_PROB_SIZE) = I(l);
		}
		for(int l=0; l < NLINES; ++l){
			x.segment(l*CL_PROB_SIZE + 1, NSEGMENTS) = Tox.row(l);
			x(l*CL_PROB_SIZE + NSEGMENTS) = -1e-6;
		}
		//
		// Set up least squares problem
		//
		int NOBS = NGUIDES*NLINES - VAL_SIZE;
		x = Barrier(AtA, Atb, x, SmoothMat, t0, MU, TOL, ALPHA, BETA, NSEGMENTS, NLINES, NGENES, NOBS, CL_PROB_SIZE, LAMBDA_G, LAMBDA_Smooth);
		//
		// Extract parameters from solution
		//
		for(int l=0; l < NLINES; ++l){
			I(l) = x(l*CL_PROB_SIZE);
		}
		//
		for(int s=0; s < NSEGMENTS; ++s){
			for(int l=0; l < NLINES; ++l){
				Tox(l, s) = x(s + l*CL_PROB_SIZE + 1);
			}
		}
		//
		for(int l=0; l < NLINES; ++l){
			G.col(l) = x.segment(NSEGMENTS + l*CL_PROB_SIZE + 1, NGENES) + x.tail(NGENES);
		}
		//
		B = x.tail(NGENES);
		//
		// Calibrate Ceff
		//
		for(int l=0; l < NLINES; ++l){
			for(int g=0; g < NGUIDES; ++g){
				for(int s=0; s < NSEGMENTS; ++s){
					if((Ceff(g, l) >= QuantileMat(l, 2*s)) && (Ceff(g, l) < QuantileMat(l, 2*s + 1))){
						double intercept = 0;
						for(int c=0; c < s; ++c){
							intercept += (QuantileMat(l, 2*c + 1) - QuantileMat(l, 2*c))*Tox(l, c);
						}
						CeffTox(g, l) = Tox(l, s)*(Ceff(g, l) - QuantileMat(l, 2*s)) + intercept;
					}
				}
			}
		}
		//
		// Report error
		//
		Dhat = Off.rowwise().replicate(Ceff.cols()) + I.colwise().replicate(Q.size()) + Q.asDiagonal()*(Mnew*G  + CeffTox);
		train_error = 0;
		val_error = 0;
		temp_error = 0;
		//
		for(int c=0; c < NREPS; ++c){
			temp_error = ((1.0/std::sqrt(reps_per_cl(ColCl(c))))*(D.col(c) - Dhat.col(ColCl(c)))).cwiseProduct(Val.col(c)).squaredNorm();
			val_error +=  temp_error / Val_Set.col(c).sum();
			train_error += (((1.0/std::sqrt(reps_per_cl(ColCl(c))))*(D.col(c) - Dhat.col(ColCl(c)))).squaredNorm() - temp_error);
		}
		val_error = val_error / NLINES;
		train_error = train_error / (NGUIDES*NLINES);
		Rcpp::Rcout << "\nTraining error is " << train_error << std::endl;
		Rcpp::Rcout << "Test error is " << val_error << std::endl;
		Rcpp::Rcout << "Smoothness is " << (Smoothness*(Tox.transpose())).squaredNorm() << std::endl;
		Rcpp::Rcout << "Gene sol size is " << G.squaredNorm() << std::endl;
		//
		if(NLINES >= 5 && fit_efficacy == 1){
			//
			// Solve least squares problem for all guides
			//
			Rcpp::Rcout << "\nSolving for guide efficacies and offsets" << std::endl;
			Eigen::MatrixXd mat = (Mnew*G + CeffTox);
			Eigen::VectorXd vec(Eigen::Map<Eigen::VectorXd>(mat.data(), mat.cols()*mat.rows()));
			val_sum = 0;
			//
			// Instantiate big sparse matrix
			//
			A = Eigen::SparseMatrix<double>(NGUIDES*NREPS - VAL_SIZE, 2*NGUIDES);
			AtripletList.clear();
			AtripletList.reserve(2*NGUIDES*NREPS);
			for(int index=0; index < NGUIDES*NREPS; ++index){
				// Set guide / cell line counter
				int gl = (index % NGUIDES) + NGUIDES*ColCl(index / NGUIDES);
				//
				if(Val_Set(index % NGUIDES, index / NGUIDES) == 0){
					// Add offset term
					AtripletList.push_back(T(val_sum, gl % NGUIDES, 1));
					// Add efficacy term
					AtripletList.push_back(T(val_sum, NGUIDES + (gl % NGUIDES), vec(gl)));
					// Increment val_sum
					val_sum += 1;
				}
			}
			A.setFromTriplets(AtripletList.begin(), AtripletList.end());
			A.makeCompressed();
			Rcpp::Rcout << "Instantiated model matrix" << std::endl;
			//
			// Reweight A
			//
			A = rep_weight.asDiagonal()*A;
			//
			// Compute A'A
			//
			AtA = A.transpose()*A;
			//
			// Set target vector
			//
			b = Eigen::Map<Eigen::VectorXd>(D.data(), D.rows()*D.cols());
			for(int c=0; c < NREPS; ++c){
				b.segment(c*NGUIDES, NGUIDES) -= I(ColCl(c))*Eigen::VectorXd::Ones(NGUIDES);
			}
			new_b = Eigen::VectorXd::Zero(NGUIDES*NREPS - VAL_SIZE);
			val_sum = 0;
			for(int gl=0; gl < NGUIDES*NREPS; ++gl){
				if(Val_Set(gl % NGUIDES, gl / NGUIDES) == 0){
					//
					new_b(val_sum) = b(gl);
					// Increment val_sum
					val_sum += 1;
				}
			}
			//
			// Reweight new_b
			//
			new_b = rep_weight.asDiagonal()*new_b;
			//
			Atb = A.transpose()*new_b;
			//
			// Instantiate guide solution vector, requiring that it is feasible, i.e. all efficacies are between 0 and 1
			//
			x = Eigen::VectorXd::Zero(2*NGUIDES);
			x.segment(0, NGUIDES) = Off;
			x.segment(NGUIDES, NGUIDES) = Q;
			//
			// Set up least squares problem
			//
			x = Barrier2(AtA, Atb, x, t02, MU, TOL, ALPHA, BETA, NGUIDES, NOBS, LAMBDA_Off);
			//
			// Extract parameters
			//
			Off = x.segment(0, NGUIDES);
			Q = x.segment(NGUIDES, NGUIDES);
		}
		//
		// Increment iter
		//
		iter += 1;
	}
	//
	Eigen::MatrixXd Val_Set = Val.cwiseProduct(Eigen::MatrixXd::Ones(NGUIDES, NREPS));
	//
	// Report error
	//
	Dhat = Off.rowwise().replicate(Ceff.cols()) + I.colwise().replicate(Q.size()) + Q.asDiagonal()*(Mnew*G  + CeffTox);
	train_error = 0;
	val_error = 0;
	temp_error = 0;
	//
	for(int c=0; c < NREPS; ++c){
		temp_error = ((1.0/std::sqrt(reps_per_cl(ColCl(c))))*(D.col(c) - Dhat.col(ColCl(c)))).cwiseProduct(Val.col(c)).squaredNorm();
		val_error +=  temp_error / Val_Set.col(c).sum();
		train_error += (((1.0/std::sqrt(reps_per_cl(ColCl(c))))*(D.col(c) - Dhat.col(ColCl(c)))).squaredNorm() - temp_error);
	}
	val_error = val_error / NLINES;
	train_error = train_error / (NGUIDES*NLINES);
	Rcpp::Rcout << "\nTraining error is " << train_error << std::endl;
	Rcpp::Rcout << "Test error is " << val_error << std::endl;
	Rcpp::Rcout << "Smoothness is " << (Smoothness*(Tox.transpose())).squaredNorm() << std::endl;
	Rcpp::Rcout << "Gene sol size is " << G.squaredNorm() << std::endl;
	//
	// Log error
	//
	loss = train_error + 0.5*LAMBDA_Smooth*(1.0 / NSEGMENTS)*(Smoothness*(Tox.transpose())).squaredNorm() + 0.5*LAMBDA_G*(1.0 / (G.rows()*G.cols()))*G.squaredNorm();
	std::string lossMessage = std::to_string(loss);
	std::string lossGMessage = "\t" + std::to_string(0.5*LAMBDA_G*(1.0 / (G.rows()*G.cols()))*G.squaredNorm());
	std::string lossSMessage = "\t" + std::to_string(0.5*LAMBDA_Smooth*(1.0 / NSEGMENTS)*(Smoothness*(Tox.transpose())).squaredNorm());
	std::string trainingErrorMessage = "\t" + std::to_string(train_error);
	std::string validationErrorMessage = "\t" + std::to_string(val_error);
	//
	std::string message = lossMessage + lossGMessage + lossSMessage + trainingErrorMessage + validationErrorMessage;
	writeTextToLog(message, logfileName);
	//
	Rcpp::Rcout << "Converting G" << std::endl;
	Rcpp::NumericVector rgenes_to_remove = Rcpp::NumericVector(genes_to_remove.sum());
	int counter = 0;
	for(int i=0; i < rG.nrow(); ++i){
		if(!genes_to_remove(i)){
			for(int j=0; j < G.cols(); ++j){
				rG(i, j) = G(counter, j);
			}
			counter += 1;
		} else {
			for(int j=0; j < G.cols(); ++j){
				rG(i, j) = NA_REAL;
			}
		}
	}
	int gene_to_remove_counter = 0;
	for(int i=0; i < genes_to_remove.size(); ++i){
		if(genes_to_remove(i) == 1){
			rgenes_to_remove[gene_to_remove_counter] = i+1;
		}
	}
	// For Tox
	Rcpp::Rcout << "Converting T" << std::endl;
	Rcpp::NumericMatrix rT = Rcpp::NumericMatrix(NLINES, NSEGMENTS);
	for(int i=0; i < Tox.rows(); ++i){
		for(int j=0; j < Tox.cols(); ++j){
			rT(i,j) = Tox(i,j);
		}
	}
	// For I
	Rcpp::Rcout << "Converting I" << std::endl;
	Rcpp::NumericVector rI = Rcpp::NumericVector(NLINES);
	for(int w=0; w < I.size(); ++w){
		rI[w] = I(w);
	}
	// For Q
	Rcpp::Rcout << "Converting Q" << std::endl;
	Rcpp::NumericVector rOff = Rcpp::NumericVector(Off.size());
	for(int w=0; w < rQ.nrow(); ++w){
		rQ(w, 2) = Q(rQ(w, 0));
	}
	for(int w=0; w < Off.size(); ++w){
		rOff[w] = Off(w);
	}
	// For B
	Rcpp::Rcout << "Converting B" << std::endl;
	Rcpp::NumericVector rB = Rcpp::NumericVector(NGENES);
	for(int w=0; w < B.size(); ++w){
		rB[w] = B(w);
	}
	Rcpp::Rcout << "Converting Dhat" << std::endl;
	Rcpp::NumericMatrix rDhat = Rcpp::NumericMatrix(NGUIDES, NLINES);
	for(int i=0; i < Dhat.rows(); ++i){
		for(int j=0; j < Dhat.cols(); ++j){
			rDhat(i, j) = Dhat(i, j);
		}
	}
	Rcpp::Rcout << "Converting CeffTox" << std::endl;
	Rcpp::NumericMatrix rCeffTox = Rcpp::NumericMatrix(NGUIDES, NLINES);
	for(int i=0; i < CeffTox.rows(); ++i){
		for(int j=0; j < CeffTox.cols(); ++j){
			rCeffTox(i, j) = CeffTox(i, j);
		}
	}
	Rcpp::Rcout << "Converting Ceff" << std::endl;
	Rcpp::NumericMatrix rCeff = Rcpp::NumericMatrix(NGUIDES, NLINES);
	for(int i=0; i < Ceff.rows(); ++i){
		for(int j=0; j < Ceff.cols(); ++j){
			rCeff(i, j) = Ceff(i, j);
		}
	}
	// For validation set
	int VAL_SIZE = Val_Set.colwise().sum().sum();
	Rcpp::NumericMatrix rVal = Rcpp::NumericMatrix(VAL_SIZE, 3);
	int row_index = 0;
	for(int i=0; i < Val.outerSize(); ++i){
		for(Eigen::SparseMatrix<double>::InnerIterator it(Val, i); it; ++it){

			if(it.value() == 1){
				rVal(row_index, 0) = it.row();
				rVal(row_index, 1) = it.col();
				rVal(row_index, 2) = 1;
				row_index += 1;
			}

		}
	}
	//
	return Rcpp::List::create(Rcpp::Named("ge_fit") = rG,
								Rcpp::Named("ce_slopes") = rT,
									Rcpp::Named("intercept_fit") = rI,
									Rcpp::Named("efficacy_fit") = rQ,
									Rcpp::Named("offset_fit") = rOff,
									Rcpp::Named("quantiles") = rQuantileMat,
									Rcpp::Named("sg_fit") = rDhat,
									Rcpp::Named("gm_fit") = rB,
									Rcpp::Named("ce_fit") = rCeffTox,
									Rcpp::Named("effective_cuts") = rCeff,
									Rcpp::Named("train_error") = train_error,
									Rcpp::Named("val_error") = val_error,
									Rcpp::Named("val_set") = rVal,
									Rcpp::Named("duplicate_genes") = rgenes_to_remove);
//
}