#!/bin/bash

##NECESSARY JOB SPECIFICATIONS                                                                                                                        
#SBATCH --partition=shared                                                                                                                            
#SBATCH --job-name=getdmu         #Set the job name to "JobExample1"                                                                                    
#SBATCH --time=0-0:05:00          #Set the wall clock limit to 1hr and 30min                                                                         
#SBATCH --ntasks=1               #Request 1 task                                                                                                      
#SBATCH --nodes=1                                                                                                                                     
#SBATCH --cpus-per-task=1                                                                                                                             
#SBATCH --mem=2000M                                                                                                                                  
#SBATCH -A tam116                                                                                                                                     
#SBATCH --output=out.%j        



python get_mu.py
