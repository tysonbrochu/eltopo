//
//  bicgstab.h
//  vortexsheet_project
//
//  Created by Tyson Brochu on 12-01-09.
//  Copyright 2012 University of British Columbia. All rights reserved.
//

#ifndef vortexsheet_project_bicgstab_h
#define vortexsheet_project_bicgstab_h

#include <blas_wrapper.h>
#include <linear_operator.h>
#include <krylov_solvers.h>

struct BiCGSTAB_Solver
{
   double tolerance_factor;
   unsigned int max_iterations;
   double residual_norm; // we use the infinity norm
   unsigned int iteration;
   KrylovSolverStatus status;
   std::vector<double> r, r_hat, s, v, p, x;
   
   BiCGSTAB_Solver(void)
   : tolerance_factor(1e-9), max_iterations(1000), residual_norm(0), iteration(0), status(KRYLOV_CONVERGED), r(0), r_hat(0), s(0), v(0), p(0), x(0)
   {}
   
   inline KrylovSolverStatus solve(const LinearOperator &A, 
                                   const double *rhs, 
                                   double *result,
                                   const LinearOperator *preconditioner=0, 
                                   bool use_given_initial_guess=false);
   
};


inline KrylovSolverStatus BiCGSTAB_Solver::solve(const LinearOperator &lin_op, 
                                                 const double *rhs, 
                                                 double *result,
                                                 const LinearOperator *preconditioner, 
                                                 bool use_given_initial_guess )
{
   const int n=lin_op.m;
   assert(lin_op.n==n);
   assert(preconditioner==0 || (preconditioner->m==n && preconditioner->n==n));
   if((int)s.size()!=n)
   {
      r.resize(n);
      r_hat.resize(n);
      s.resize(n);
   }
   
   // convergence tolerance
   double tol=tolerance_factor*BLAS::abs_max(n, rhs);
  
   if ( !use_given_initial_guess )
   {
      BLAS::set_zero(n, result);
   }
   
   // r0 = b - A*x0
   lin_op.apply_and_subtract( result, rhs, &r[0] );

   residual_norm=BLAS::abs_max(r);
   if(residual_norm==0) return status=KRYLOV_CONVERGED;

   // r_hat0 = r0
   r_hat = r;
   
   // rho = alpha = omega = 1
   double rho = 1.0;
   double alpha = 1.0;
   double omega = 1.0;
   
   // v = p = 0
   v.resize(n,0);
   p.resize(n,0);
   
   std::vector<double> product;
   std::vector<double> t;
   for(iteration=1; iteration<max_iterations; ++iteration)
   {
      double rho_prev = rho;
      rho = BLAS::dot(r_hat,r);
            
      double beta = (rho/rho_prev) * (alpha/omega);
      
      if( rho!=rho || beta!=beta )
      {
         std::cout << "BiCGSTAB_Solver: NaN detected, restarting" << std::endl;
         BLAS::set_zero(n, result);
         lin_op.apply_and_subtract( result, rhs, &r[0] );
         iteration = 1;
         rho = 1.0;
         alpha = 1.0;
         omega = 1.0;
         continue;
      }
      
      // p = r + beta*(p - omega*v);
      BLAS::add_scaled( -omega, v, p );
      product = r;
      BLAS::add_scaled( beta, p, product );
      p = product;
      
      // v = A*p
      lin_op.apply( p, v );

      // alpha = rho / (r_hat,v)
      alpha = rho / BLAS::dot(r_hat,v);
      
      // s = r - alpha * v
      s = r;
      BLAS::add_scaled( -alpha, v, s );

      // t = A*s
      lin_op.apply( s, t );
      
      // omega = (t,s)/(t,t)
      omega = BLAS::dot(t,s) / BLAS::dot(t,t);
      
      // x = x + alpha*p + omega*s
      BLAS::add_scaled( n, alpha, &p[0], result );
      BLAS::add_scaled( n, omega, &s[0], result );
      
      // check convergence
      residual_norm=BLAS::abs_max(r);
      
      if(residual_norm<=tol) 
      {
         return status=KRYLOV_CONVERGED;
      }
      
      // r = s - omega * t
      BLAS::add_scaled( -omega, t, s );
      r = s;
   }
   
   return status=KRYLOV_EXCEEDED_MAX_ITERATIONS;
   
}


#endif


