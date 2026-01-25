/*
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <fstream>
#include <iterator>
#include <iostream>
#include <unordered_set>
#include <vector>
#include <numeric>
#include <limits>
#include <cstdlib>
#include <boost/math/distributions.hpp>

#include <savvy/reader.hpp>
#include <savvy/writer.hpp>

#include "tcdf.hpp"

template<typename T>
T square(T v)
{
  return v * v;
}

float strtof(const std::string& s)
{
  return s.empty() || s[0] == 'N' ? std::numeric_limits<float>::quiet_NaN() : std::atof(s.c_str());
}

std::vector<std::string> str_split(std::string s, const std::string& delimiter) {
  std::vector<std::string> ret;
  size_t pos = 0;
  while ((pos = s.find(delimiter)) != std::string::npos)
  {
    ret.push_back(s.substr(0, pos));
    s.erase(0, pos + delimiter.length());
  }
  ret.push_back(s);

  return ret;
}

std::unordered_set<std::string> split_file_to_set(const char* in)
{
  std::unordered_set<std::string> ret;

  std::string s;
  std::ifstream ifs(in);
  while (std::getline(ifs, s))
  {
    ret.insert(s);
  }

  return ret;
}

std::vector<std::string> split_file_to_vector(const char* in, std::size_t size_hint = 100)
{

  std::vector<std::string> ret;
  ret.reserve(size_hint);

  std::string s;
  std::ifstream ifs(in);
  while (std::getline(ifs, s))
  {
    ret.emplace_back(s);
  }

  return ret;
}

int main(int argc, char** argv)
{
  auto chip_genes = split_file_to_set("data/whitelist_filter_20230531/NEJM_2017_genes_01262020_first_col.MLL_fix.txt");
  auto missense_variants = split_file_to_set("data/whitelist_filter_20230531/CHIP_missense_vars_cv_04102022.txt");
  auto splice_genes = split_file_to_set("data/whitelist_filter_20230531/CHIP_splice_vars_agb_01262020.txt");
  auto lof_genes = split_file_to_set("data/whitelist_filter_20230531/CHIP_nonsense_FS_vars_agb_01262020.txt");

  std::vector<float> ages;
  {
    std::vector<std::string> ages_str;
    std::ifstream is("data/merged_ad.nwd_ids.age.txt");
    bool b = is.good();
    b = is.good();
    std::istream_iterator<std::string> s(is), e;
    for ( ; s != e; ++s)
      ages.push_back(strtof(*s));
  }

  savvy::reader input_file(argv[1]);

  if (input_file.samples().size() != ages.size())
    return std::cerr << "Error: age file should contain lines for each sample in the VCF and with the same order\n", EXIT_FAILURE;
  //input_file.reset_bounds(savvy::genomic_region("chr2", 25234372, 25234373)); R882
  //input_file.reset_bounds(savvy::genomic_region("chr2", 25247132, 25247133));
  //input_file.reset_bounds(savvy::genomic_region("chr21",43092956, 43107570));
  input_file.reset_bounds(savvy::genomic_region("chr17", 7669662, 7669663));
  bool b = input_file.good();

  auto hdrs = input_file.headers();
  hdrs.emplace_back("FILTER", "<ID=OFF_TARGET, Description=\"Not an exonic of splicing variant in target gene\">");
  hdrs.emplace_back("FILTER", "<ID=GERM, Description=\"Failed germline test\">");
  hdrs.emplace_back("FILTER", "<ID=AGE_ASSOC, Description=\"Failed age association test\">");
  hdrs.emplace_back("INFO", "<ID=NC,Number=1,Type=Integer,Description=\"Number of CHIP carriers\">");
  hdrs.emplace_back("INFO", "<ID=NCA,Number=1,Type=Integer,Description=\"Number of CHIP carriers with known age\">");
  hdrs.emplace_back("INFO", "<ID=AGE_P,Number=1,Type=Float,Description=\"Age association test log10 p-value\">");
  hdrs.emplace_back("INFO", "<ID=GERM_P,Number=1,Type=Float,Description=\"Germline test p-value\">");
  hdrs.emplace_back("INFO", "<ID=KNOWN_MIS,Number=0,Type=Flag,Description=\"Variant is a known missense CHIP mutation\">");
  hdrs.emplace_back("FORMAT","<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fractions with samples having less than two supporting (i.e., ALT) reads set to zero\">");

  savvy::writer output_file(argv[2], savvy::file::format::vcf, hdrs, input_file.samples());

  std::vector<float> vaf(input_file.samples().size());
  savvy::variant rec;
  while (output_file && input_file >> rec)
  {
    auto flt = rec.filters();
    std::string s;
    rec.get_info("Gene.refGene", s);
    auto genes = str_split(s, "\\x3b");
    rec.get_info("Func.refGene", s);
    auto funcs = str_split(s, "\\x3b");

//    if (genes.size() != funcs.size())
//    {
//      std::cerr << "Notice: gene and func vector lengths do not match: " << rec.chrom() <<":" << rec.pos() << "-" << rec.pos() << std::endl;
//    }

    std::size_t idx = 0;
    for ( ; idx < genes.size(); ++idx)
    {
      std::size_t func_idx = funcs.size() == genes.size() ? idx : 0;
      if (chip_genes.find(genes[idx]) != chip_genes.end() && (funcs[func_idx] == "exonic" | funcs[func_idx] == "splicing"))
        break;
    }

    bool known_mis = false;
    if (idx >= genes.size())
    {
      // TODO: mark fail
      if (flt.size() == 1 && flt[0] == "PASS")
        flt.clear();
      flt.push_back("OFF_TARGET");
      //std::cerr << rec.chrom() << ":" << rec.pos() << "_" << rec.ref() << "/" << (rec.alts().empty() ? "" : rec.alts()[0]) << "\tNot in gene list" << std::endl;

    }
    else
    {
      rec.get_info("ExonicFunc.refGene", s);
      if (s == "nonsynonymous_SNV")
      {
        rec.get_info("AAChange.refGene", s);

        auto aa_change = str_split(s, ",");
        for (std::size_t i = 0; i < aa_change.size(); ++i)
        {
          auto aa = str_split(aa_change[i], ":");
          if (aa.size() != 5 || aa[4].size() < 3)
          {
            std::cerr << "Error: cannot parse AAChange.refGene: " << rec.chrom() << ":" << rec.pos() << std::endl;
          }
          else
          {
            if (missense_variants.find(aa[0] + "\t" + aa[4].substr(2)) != missense_variants.end())
            {
              known_mis = true;
              break;
            }
          }
        }
      }
    }

    std::vector<std::int16_t> ad;
    rec.get_format("AD", ad);

    if (ad.size() != vaf.size() * 2)
      return std::cerr << "Error: AD has wrong length\n", EXIT_FAILURE;

    double mu = 0.;
    std::array<double, 2> age_mu{};// {non-carrier, carrier}
    std::array<std::int32_t, 2> age_n{};
    std::int32_t n_carriers = 0;
    for (std::size_t i = 0; i < vaf.size(); ++i)
    {
      float dp = float(ad[i * 2] + ad[i * 2 + 1]);
      if (std::isfinite(dp) && ad[i * 2 + 1] > 1)
      {
        vaf[i] = float(ad[i * 2 + 1]) / dp;
        mu += vaf[i];

        if (std::isfinite(ages[i]))
        {
          age_mu[1] += ages[i];
          ++age_n[1];
        }

        ++n_carriers;
      }
      else
      {
        vaf[i] = 0.;
        if (std::isfinite(ages[i]))
        {
          age_mu[0] += ages[i];
          ++age_n[0];
        }
      }
    }

    rec.set_info("NC", std::int32_t(n_carriers));
    rec.set_info("NCA", std::int32_t(age_n[1]));
    rec.set_format("VAF", vaf);

    mu = mu / n_carriers;
    age_mu[0] = age_mu[0] / age_n[0];
    age_mu[1] = age_mu[1] / age_n[1];

    if (n_carriers > 1)
    {
      double sd = 0.;
      for (std::size_t i = 0; i < vaf.size(); ++i)
      {
        if (vaf[i] > 0.f)
        {
          double e = vaf[i] - mu;
          sd += e * e;
        }
      }

      sd = std::sqrt(sd / (n_carriers - 1.));
      double t = (0.5 - mu) / (sd / std::sqrt(n_carriers));
      boost::math::students_t_distribution<double> dist(n_carriers - 1);
      double germ_pval =  boost::math::cdf(complement(dist, std::isnan(t) ? 0. : t));

      rec.set_info("GERM_P", float(std::log10(germ_pval)));
      if (germ_pval >= 0.05)
      {
        if (flt.size() == 1 && flt[0] == "PASS")
          flt.clear();
        flt.push_back("GERM");
      }
    }

    if (age_n[0] > 1 && age_n[1] > 1)
    {
      std::array<double, 2> sv{};

      for (std::size_t i = 0; i < ages.size(); ++i)
      {
        std::size_t j = vaf[i] > 0.f ? 1 : 0;
        if (std::isfinite(ages[i]))
        {
          double e = ages[i] - age_mu[j];
          sv[j] += e * e; //std::max(e * e, std::numeric_limits<double>::min());
        }
      }

      for (std::size_t j = 0; j < 2; ++j)
        sv[j] = sv[j] / (age_n[j] - 1.);

      double t = (age_mu[1] - age_mu[0]) / std::sqrt(sv[0] / age_n[0] + sv[1] / age_n[1]);
      double dof = square((sv[0] / age_n[0] + sv[1] / age_n[1])) / (square(sv[0] / age_n[0]) / (age_n[0] - 1.) + square(sv[1] / age_n[1]) / (age_n[1] - 1.));
      boost::math::students_t_distribution<double> dist(dof);
      double age_pval =  boost::math::cdf(complement(dist, std::isnan(t) ? 0. : t));

      rec.set_info("AGE_P", float(std::log10(age_pval)));
      if (age_pval >= 0.05)
      {
        if (flt.size() == 1 && flt[0] == "PASS")
          flt.clear();
        flt.push_back("AGE_ASSOC");
      }
    }

    if (known_mis)
    {
      rec.set_info("KNOWN_MIS", std::vector<std::int8_t>());
      flt = {"PASS"};
    }

    int ns;
    rec.get_info("NS", ns);
    for (std::string info : {"DP","ECNT","NLOD","N_ART_LOD","POP_AF","P_CONTAM","P_GERMLINE","TLOD"})
    {
      float v;
      rec.get_info(info, v);
      rec.set_info(info, v / ns); // These Mutect2 INFO fields were summed during `bcftools merge`. Make them averages.
    }

    rec = savvy::site_info(rec.chrom(), rec.pos(), rec.ref(), rec.alts(), rec.id(), rec.qual(), flt, rec.info_fields());

    output_file << rec;
    //std::cerr << s << "\n";
  }

  return 0;
}
