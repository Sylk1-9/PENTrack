#include <cmath>
#include <fstream>
#include <iostream>

#include "boost/format.hpp"
#include <boost/iterator/zip_iterator.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "boost/multi_array.hpp"

#include "globals.h"
#include "ecurrentFields.h"


TECurrentField::TECurrentField(const std::string sft, const std::string &_It) {

  // TCustomBField::TCustomBField(const std::string &_Bx, const std::string &_By, const std::string &_Bz){
  
  tvar = std::unique_ptr<double>(new double(0.0));
  exprtk::symbol_table<double> symbol_table;
  symbol_table.add_variable("t",*tvar);
  symbol_table.add_constants();
  exprtk::parser<double> parser;
  std::string expr{_It,};
  Itexpr.register_symbol_table(symbol_table);
  if (not parser.compile(expr, Itexpr)){
    throw std::runtime_error(exprtk::parser_error::to_str(parser.get_error(0).mode) + " while parsing Custom It Current formula '" + expr + "': " + parser.get_error(0).diagnostic);
  }

  
  boost::filesystem::path ft(sft);
  ft = boost::filesystem::absolute(ft, configpath.parent_path());

  std::string line;
  std::vector<std::string> line_parts;
  WireSegment segment;

  std::ifstream FINstream(ft.string(), std::ifstream::in);
  boost::iostreams::filtering_istream FIN;
  if (boost::filesystem::extension(ft) == ".bz2"){
    FIN.push(boost::iostreams::bzip2_decompressor());
  }
  else if (boost::filesystem::extension(ft) == ".gz"){
    FIN.push(boost::iostreams::gzip_decompressor());
  }
  FIN.push(FINstream);
  if (!FINstream.is_open() or !FIN.is_complete()){
    throw std::runtime_error("Could not open " + ft.string());
  }
  std::cout << "\nReading " << ft << "\n";

  // Read in file data
  int lineNum = 0;
  while (getline(FIN,line)){
    lineNum++;
    if (line.substr(0,1) == "%" || line.substr(0,1) == "#") continue;     // Skip commented lines
    boost::split(line_parts, line, boost::is_any_of("\t, "), boost::token_compress_on); //Delineate tab, space, commas
    
    if (line_parts.size() < 3){
      throw std::runtime_error((boost::format("Error reading line %1% of file %2%") % lineNum % ft.string()).str());
    }
    
    segment.x1 = std::stod(line_parts[0], nullptr);
    segment.y1 = std::stod(line_parts[1], nullptr);
    segment.z1 = std::stod(line_parts[2], nullptr);
    segment.x2 = std::stod(line_parts[3], nullptr);
    segment.y2 = std::stod(line_parts[4], nullptr);
    segment.z2 = std::stod(line_parts[5], nullptr);
    
    // segment.current = current;
    wireSegments.push_back(segment);    
  }

  if (wireSegments.empty()) {
    throw std::runtime_error("No data read from " + ft.string());
  }
  
}



void TECurrentField::BField(const double x, const double y, const double z, const double t, double B[3], double dBidxj[3][3]) const {
  *tvar = t;
  double current =  Itexpr.value();

  double bx = 0.0;
  double by = 0.0;
  double bz = 0.0;

  double dbx_dxi[3] = {0.0};
  double dby_dxi[3] = {0.0};
  double dbz_dxi[3] = {0.0};

  for (const auto& segment : wireSegments) {
    double x1 = segment.x1;
    double y1 = segment.y1;
    double z1 = segment.z1;
    double x2 = segment.x2;
    double y2 = segment.y2;
    double z2 = segment.z2;
    // double current = segment.current;

    double dx = x2 - x1;
    double dy = y2 - y1;
    double dz = z2 - z1;

    double r = std::sqrt(dx * dx + dy * dy + dz * dz);
    double r3 = r * r * r;

    double dlx = mu0/(4*pi) * current * (dy * (z - z1) - dz * (y - y1)) / r3;
    double dly = mu0/(4*pi) * current * (dz * (x - x1) - dx * (z - z1)) / r3;
    double dlz = mu0/(4*pi) * current * (dx * (y - y1) - dy * (x - x1)) / r3;

    bx += dlx;
    by += dly;
    bz += dlz;

    // First derivatives
    if (dBidxj != nullptr){
      double dr_dxi[3] = {0.0};
      dr_dxi[0] = (x - x1) / r;
      dr_dxi[1] = (y - y1) / r;
      dr_dxi[2] = (z - z1) / r;

      for (int i = 0; i < 3; i++) {
	dbx_dxi[i] += mu0/(4*pi) * current * (dy * dr_dxi[i] - dz * dr_dxi[(i + 1) % 3]) / r3;
	dby_dxi[i] += mu0/(4*pi) * current * (dz * dr_dxi[i] - dx * dr_dxi[(i + 1) % 3]) / r3;
	dbz_dxi[i] += mu0/(4*pi) * current * (dx * dr_dxi[i] - dy * dr_dxi[(i + 1) % 3]) / r3;
      }
    }
  }

  B[0] = bx;
  B[1] = by;
  B[2] = bz;

  if (dBidxj != nullptr){

    dBidxj[0][0] = dbx_dxi[0];
    dBidxj[0][1] = dbx_dxi[1];
    dBidxj[0][2] = dbx_dxi[2];

    dBidxj[1][0] = dby_dxi[0];
    dBidxj[1][1] = dby_dxi[1];
    dBidxj[1][2] = dby_dxi[2];

    dBidxj[2][0] = dbz_dxi[0];
    dBidxj[2][1] = dbz_dxi[1];
    dBidxj[2][2] = dbz_dxi[2];
  }
}
