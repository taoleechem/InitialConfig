/*
 * =====================================================================================
 *
 *       Filename:  configs.h
 *
 *    Description:  Generate some configurations of 2 molecules for the automatically
 *                  search of chemical reaction path
 *
 *        Version:  1.0
 *        Created:  2016年03月21日 14时39分14秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Tao Li (), taoleechem@outlook.com
 *   Organization:  Nanjing University
 *
 * =====================================================================================
 */

#ifndef _CONFIGS_H_
#define _CONFIGS_H_
#include <iostream>
#include <vector>
void Do_GenerateFunction_Program_FromFile(std::string filename);
void Do_RandomGenerate_FromFile(std::string filename);
void WriteWaters(std::string filename, std::string tip5p_file, int filename_single_num);

#endif
