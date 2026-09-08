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
#include <regex>
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

bool write_detail(std::ostream& os, savvy::site_info& site, const std::string& gene_detail, const std::string& aa_change, float germ_p, float age_p, std::int32_t nc, std::int32_t nca)
{
  os << site.chrom()
    << "\t" << site.pos()
    << "\t" << site.ref()
    << "\t" << (site.alts().empty() ? "." : site.alts()[0])
    << "\t" << (gene_detail.empty() ? "." : gene_detail)
    << "\t" << (aa_change.empty() ? "." : aa_change)
    << "\t" << germ_p
    << "\t" << age_p
    << "\t" << nc
    << "\t" << nca;
  return os.put('\n').good();
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
  const int min_alt_depth = 3;
  const int min_alt_depth_u2af1 = 5;
  const float min_vaf = 0.02f;
  const bool apply_germ_and_age_filters = false;

  std::regex lof_regex("fs|X|\\*"); // fs is probably redundant as the stop codon "X" or "*" should follow
  /*bool m2 = std::regex_search("p.X123R", lof_regex);
  m2 = std::regex_search("p.R123X", lof_regex);
  m2 = std::regex_search("p.R123*", lof_regex);
  m2 = std::regex_search("p.R123H", lof_regex);
  m2 = std::regex_search("p.R123fs*10", lof_regex);
  m2 = std::regex_search("p.R97fsX121", lof_regex);*/

#if 0
  auto chip_genes = split_file_to_set("data/whitelist_filter_20230531/NEJM_2017_genes_01262020_first_col.MLL_fix.txt");
  auto chip_gene_accessions = split_file_to_set("data/whitelist_filter_20230531/NEJM_2017_genes_01262020_nocr_mll_fix.accessions.tsv");
  //auto missense_variants = split_file_to_set("data/whitelist_filter_20230531/CHIP_missense_vars_cv_04102022.txt");
  auto missense_variants = split_file_to_set("data/whitelist_filter_20230531/NEJM_2017_genes_01262020.MLL_fix.flatten.cleaned.additional_missense_v3.tsv");
  //auto lof_genes = split_file_to_set("data/whitelist_filter_20230531/CHIP_nonsense_FS_vars_agb_01262020.MLL_fix.txt");
  auto lof_genes = split_file_to_set("data/whitelist_filter_20230531/NEJM_2017_genes_01262020_nocr_mll_fix.lof.tsv");
  //auto splice_genes = split_file_to_set("data/whitelist_filter_20230531/CHIP_splice_vars_agb_01262020.txt");
  auto splice_genes_vec = split_file_to_vector("data/whitelist_filter_20230531/NEJM_2017_genes_01262020_nocr_mll_fix.splice.tsv");
  std::unordered_set<std::string> splice_genes;
  splice_genes.reserve(splice_genes_vec.size());
  for (auto it = splice_genes_vec.begin(); it != splice_genes_vec.end(); ++it)
  {
    auto fields = str_split(*it, "\t");
    if (fields.size() != 2)
      return std::cerr << "Error: could not parse splicing genee list\n", EXIT_FAILURE;
    splice_genes.insert(fields[1]);
  }
#endif
  std::vector<float> ages;
  {
    std::vector<std::string> ages_str;
    std::ifstream is("data/merged_ad.nwd_ids.age.txt");
    //std::ifstream is("data/filter_param.txt");
    bool b = is.good();
    b = is.good();
    std::istream_iterator<std::string> s(is), e;
    for ( ; s != e; ++s)
      ages.push_back(strtof(*s));
  }

  savvy::reader input_file(argv[1]);

  if (input_file.samples().size() != ages.size())
    return std::cerr << "Error: age file should contain lines for each sample in the VCF and with the same order\n", EXIT_FAILURE;
  //input_file.reset_bounds(savvy::genomic_region("chr2", 25234372, 25234373)); //R882
  //input_file.reset_bounds(savvy::genomic_region("chr2", 25247132, 25247133));
  //input_file.reset_bounds(savvy::genomic_region("chr21",43092956, 43107570));
  //input_file.reset_bounds(savvy::genomic_region("chr17", 7669662, 7669663));
  //input_file.reset_bounds(savvy::genomic_region("chr12", 49019423, 49060794));
  //input_file.reset_bounds(savvy::genomic_region("chr12", 49051672, 49060794));
  //input_file.reset_bounds(savvy::genomic_region("chr20", 32433283, 32433447));
  //input_file.reset_bounds(savvy::genomic_region("chr2", 25234418, 25235710));
  bool b = input_file.good();

  auto hdrs = input_file.headers();
  //hdrs.emplace_back("FILTER", "<ID=OFF_TARGET, Description=\"Not an exonic of splicing variant in target gene\">");
  if (apply_germ_and_age_filters)
  {
    hdrs.emplace_back("FILTER", "<ID=GERM, Description=\"Failed germline test\">");
    hdrs.emplace_back("FILTER", "<ID=AGE_ASSOC, Description=\"Failed age association test\">");
  }
  hdrs.emplace_back("INFO", "<ID=NC,Number=1,Type=Integer,Description=\"Number of CHIP carriers\">");
  hdrs.emplace_back("INFO", "<ID=NCA,Number=1,Type=Integer,Description=\"Number of CHIP carriers with known age\">");
  hdrs.emplace_back("INFO", "<ID=AGE_P,Number=1,Type=Float,Description=\"Age association test log10 p-value\">");
  hdrs.emplace_back("INFO", "<ID=GERM_P,Number=1,Type=Float,Description=\"Germline test p-value\">");
  hdrs.emplace_back("INFO", "<ID=AGE_M_DIFF,Number=1,Type=Float,Description=\"Difference between mean age of carriers and non-carriers\">");
  hdrs.emplace_back("INFO", "<ID=VAF_M,Number=1,Type=Float,Description=\"Mean VAF\">");
//  hdrs.emplace_back("INFO", "<ID=KNOWN_CHIP,Number=0,Type=Flag,Description=\"Variant is a known missense or LOF CHIP mutation\">");
//  hdrs.emplace_back("INFO", "<ID=LOF,Number=0,Type=Flag,Description=\"Variant is nonsense, nonstop, or frameshif and in loss of function gene list\">");
//  hdrs.emplace_back("INFO", "<ID=SPLICE,Number=0,Type=Flag,Description=\"Is a predicted splicing variant and in splicing gene list\">");
  hdrs.emplace_back("FORMAT","<ID=VAF,Number=1,Type=Float,Description=\"Variant allele fractions with zero thresholding applied for min ALT depth and min VAF\">");

  savvy::writer output_file(argv[2], savvy::file::format::bcf, hdrs, input_file.samples());
//  shrinkwrap::bgzf::ostream detail_file(argv[3]);
//  detail_file << "CHROM\tPOS\tREF\tALT\tGENE_DETAIL\tAA_CHANGE\tGERM_P\tAGE_P\tNC\tNCA" << std::endl;
//  if (!detail_file)
//    return std::cerr << "Error: could not open detail file\n", EXIT_FAILURE;

  std::vector<float> vaf(input_file.samples().size());
  savvy::variant rec;
  while (output_file && input_file >> rec)
  {
//    if (genes.size() != funcs.size())
//    {
//      std::cerr << "Notice: gene and func vector lengths do not match: " << rec.chrom() <<":" << rec.pos() << "-" << rec.pos() << std::endl;
//    }

    std::vector<std::int16_t> ad;
    rec.get_format("AD", ad);

    if (ad.size() != vaf.size() * 2)
      return std::cerr << "Error: AD has wrong length\n", EXIT_FAILURE;

    std::string gene;
    rec.get_info("Gene.refGene", gene);
    bool is_u2af1 =  gene == "U2AF1\\x3bU2AF1L5";

    double mu = 0.;
    std::array<double, 2> age_mu{};// {non-carrier, carrier}
    std::array<std::int32_t, 2> age_n{};
    std::int32_t n_carriers = 0;
    for (std::size_t i = 0; i < vaf.size(); ++i)
    {
      float dp = float(ad[i * 2] + ad[i * 2 + 1]);
      float alt_dp = ad[i * 2 + 1];
      float v = alt_dp / dp;
      if (std::isfinite(dp) && alt_dp >= (is_u2af1 ? min_alt_depth_u2af1 : min_alt_depth) && v >= min_vaf)
      {
        vaf[i] = v;
        mu += v;

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

    if (n_carriers > 0)
      rec.set_info("VAF_M", float(mu));

    if (age_n[0] > 0 && age_n[1] > 0)
      rec.set_info("AGE_M_DIFF", float(age_mu[1] - age_mu[0]));

    double germ_pval = std::numeric_limits<double>::quiet_NaN();
    double age_pval = std::numeric_limits<double>::quiet_NaN();

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
      germ_pval =  boost::math::cdf(complement(dist, std::isnan(t) ? 0. : t));
      rec.set_info("GERM_P", float(std::log10(germ_pval)));
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
      age_pval =  boost::math::cdf(complement(dist, std::isnan(t) ? 0. : t));
      rec.set_info("AGE_P", float(std::log10(age_pval)));
    }

#if 0
    bool off_target = false;
    bool known_mis = false;
    bool is_lof = false;
    bool is_splice = false;
    bool has_splicing_func = false;
    bool has_exonic_func = false;

    std::string s;
    rec.get_info("Gene.refGene", s);
    auto genes = str_split(s, "\\x3b");
    rec.get_info("Func.refGene", s);
    auto funcs = str_split(s, "\\x3b");
    rec.get_info("ExonicFunc.refGene", s);
    auto exon_funcs = str_split(s, "\\x3b");
    std::size_t idx = 0;
    for ( ; idx < genes.size(); ++idx)
    {
      std::size_t func_idx = funcs.size() == genes.size() ? idx : 0;
      std::size_t exon_func_idx = exon_funcs.size() == genes.size() ? idx : 0;
      if (chip_genes.find(genes[idx]) != chip_genes.end())
      {
        if (funcs[func_idx] == "exonic" && exon_funcs[exon_func_idx] != "synonymous_SNV")
        {
          has_exonic_func = true;
          break;
        }
        else if (funcs[func_idx] == "splicing")
        {
          has_splicing_func = true;
          break;
        }
      }
    }

    if (idx >= genes.size())
    {
      off_target = true;
      //std::cerr << rec.chrom() << ":" << rec.pos() << "_" << rec.ref() << "/" << (rec.alts().empty() ? "" : rec.alts()[0]) << "\tNot in gene list" << std::endl;
    }
    else
    {
      std::string aa_out;
      std::string gene_detail_out;

      if (has_exonic_func)
      {
        rec.get_info("AAChange.refGene", s); // rec.get_info("ExonicFunc.refGene", s);
        if (!s.empty() && s != ".")          // s == "nonsynonymous_SNV")
        {
          auto aa_change = str_split(s, ",");
          for (std::size_t i = 0; i < aa_change.size(); ++i)
          {
            auto aa = str_split(aa_change[i], ":");
            if (aa.size() != 5 || aa[4].size() < 3)
            {
              // std::cerr << "Error: cannot parse AAChange.refGene: " << rec.chrom() << ":" << rec.pos() << std::endl;
            }
            else
            {
              if (chip_gene_accessions.find(aa[1]) != chip_gene_accessions.end())
              {
                if (aa_out.size() > 0)
                {
                  std::cerr << "Warning: duplicagted accession record - aa change\n";
                }
                else
                {
                  aa_out = aa_change[i];
                }
              }

              // if (missense_variants.find(aa[0] + "\t" + aa[4].substr(2)) != missense_variants.end())
              if (missense_variants.find(aa[0] + "\t" + aa[1] + "\t" + aa[4]) != missense_variants.end())
              {
                known_mis = true;
                // break;
              }
              else if (lof_genes.find(aa[0] + "\t" + aa[1]) != lof_genes.end() && std::regex_search(aa[4], lof_regex))
              {
                is_lof = true;
              }
            }
          }
        }
      }

      if (has_splicing_func)
      {
        rec.get_info("GeneDetail.refGene", s); // rec.get_info("ExonicFunc.refGene", s);
        if (!s.empty() && s != ".")          // s == "nonsynonymous_SNV")
        {
          auto gene_detail = str_split(s, "\\x3b");
          for (std::size_t i = 0; i < gene_detail.size(); ++i)
          {
            auto d = str_split(gene_detail[i], ":");
            if (d.size() < 2 || d[1].size() < 4 || d[1].compare(0, 4, "exon") != 0)
            {
              // std::cerr << "Error: cannot parse GeneDetail.refGene: " << rec.chrom() << ":" << rec.pos() << std::endl;
            }
            else
            {
              if (chip_gene_accessions.find(d[0]) != chip_gene_accessions.end())
              {
                if (gene_detail_out.size() > 0)
                {
                  std::cerr << "Warning: duplicagted accession record - gene detail\n";
                }
                else
                {
                  gene_detail_out = gene_detail[i];
                }
              }

              // if (missense_variants.find(aa[0] + "\t" + aa[4].substr(2)) != missense_variants.end())
              if (splice_genes.find(d[0]) != splice_genes.end())
              {
                is_splice = true;
                // break;
              }
            }
          }
        }
      }

      if (aa_out.size() + gene_detail_out.size() > 0)
      {
        if (!write_detail(detail_file, rec, gene_detail_out, aa_out, germ_pval, age_pval, n_carriers, age_n[1]))
          return std::cerr << "Error: failed to write to detail file\n", EXIT_FAILURE;
      }
    }

    if (!known_mis)
    {

    }

    /*if (n_carriers > 1 && age_n[0] > 1 && age_n[1] > 1) // Both t-tests will be applied
    {
      flt = {"PASS"};
    }

    if (n_carriers > 1 && germ_pval >= 0.05)
    {
      if (flt.size() == 1 && flt[0] == "PASS")
        flt.clear();
      flt.push_back("GERM");
    }

    if (age_n[0] > 1 && age_n[1] > 1 && age_pval >= 0.05)
    {
      if (flt.size() == 1 && flt[0] == "PASS")
        flt.clear();
      flt.push_back("AGE_ASSOC");
    }

    if (off_target)
    {
      if (flt.size() == 1 && flt[0] == "PASS")
        flt.clear();
      flt.push_back("OFF_TARGET");
    }
    else if (known_mis)
    {
      rec.set_info("KNOWN_CHIP", std::vector<std::int8_t>());
      flt = {"PASS"};
    }*/

    if (known_mis)
    {
      rec.set_info("KNOWN_CHIP", std::vector<std::int8_t>());
    }

    if (is_lof)
    {
      rec.set_info("LOF", std::vector<std::int8_t>());
    }

    if (is_splice)
    {
      rec.set_info("SPLICE", std::vector<std::int8_t>());
    }

    std::vector<std::string> flt; //= rec.filters();

    if (known_mis)
    {
      flt = {"PASS"};
    }
    else
    {
      if (!is_lof && !is_splice)
        flt.push_back("OFF_TARGET");

      if (apply_germ_and_age_filters)
      {
        if (n_carriers > 1 && germ_pval >= 0.05)
          flt.push_back("GERM");

        if (age_n[0] > 1 && age_n[1] > 1 && age_pval >= 0.05)
          flt.push_back("AGE_ASSOC");
      }

      if (flt.empty())
        flt = {"PASS"};
    }


//    int ns;
//    rec.get_info("NS", ns);
//    for (std::string info : {"DP","ECNT","NLOD","N_ART_LOD","POP_AF","P_CONTAM","P_GERMLINE","TLOD"})
//    {
//      float v;
//      rec.get_info(info, v);
//      rec.set_info(info, v / ns); // These Mutect2 INFO fields were summed during `bcftools merge`. Make them averages.
//    }
#endif

    auto flt = rec.filters();
    rec = savvy::site_info(rec.chrom(), rec.pos(), rec.ref(), rec.alts(), rec.id(), rec.qual(), flt, rec.info_fields());

    output_file << rec;
    //std::cerr << s << "\n";
  }

  return 0;
}
