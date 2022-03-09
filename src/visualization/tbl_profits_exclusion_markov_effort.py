#!/usr/bin/env python3

from profitability import *
import os

z=1
n_main=0
n_fork=0
kr=1
kl=10
N=8
Prslt="Pwin"

p_B=0
p_A=None
p_V=0
p_I=0
p_E=0
p_m=0

if "P_V" in os.environ:
    p_V = float(str(os.environ["P_V"]))
else:
    p_V=0

if "Z" in os.environ:
    z = int(str(os.environ["Z"]))
else:
    z=1

if "N_FORK" in os.environ:
    n_fork = int(os.environ["N_FORK"])
else:
    n_fork = 0

if "N_MAIN" in os.environ:
    n_main = int(os.environ["N_MAIN"])
else:
    n_main = 0

if "KR" in os.environ:
    kr = int(str(os.environ["KR"]))
else:
    kr = 1

if "KL" in os.environ:
    kl = int(str(os.environ["KL"]))
else:
    kl = 10

sympy_expr=get_sympy_expr_for_e_min_markov_blocks_effort()[0]

table = """
\\renewcommand{\\arraystretch}{1.3} % increase space between lines
\\begin{table*}[!htpb]
\\centering
\\scriptsize
\\resizebox{1.0\\linewidth}{!}{   % 0.7 for 2 colum view
\\begin{tabular}{c|c|c|c|c|c|c|c|c|}
"""

O_m = [0.05, 0.10, 0.20, 0.30, 0.33, 0.382, 0.40, ]
O_a = [0.0, 0.05, 0.10, 0.20, 0.30, 0.33, 0.382, 0.40, ]

for i in O_m:
    table += "\t& $ p_m = $ {} ".format(str(i))
#print("\\ \\hline")

for p_E in O_a:
    table += """\\\\ \\hline
    """
    table += "$ p_\\mathcal{{ E }} = {:.2f} $".format(p_E)
    for p_m in O_m:
        if 1 - (p_E + p_m) < 0:
            continue

        min_cost = e_min_markov_sympy(p_B=p_B,
                                p_E=p_E,
                                p_V=p_V,
                                p_I=p_I,
                                p_m=p_m,
                                z=z,
                                n_main=n_main,
                                n_fork=n_fork,
                                kr=kr,
                                kl=kl,
                                N=N,
                                expr=sympy_expr)

        proift_main_with_m = EV_main_markov_blocks(p_B=p_B,
                                                  p_E=p_E,
                                                  p_V=p_V,
                                                  p_I=p_I,
                                                  p_m=p_m,
                                                  z=z,
                                                  n_fork=n_fork,
                                                  n_main=n_main,
                                                  epsilon=min_cost,
                                                  kr=kr,
                                                  kl=kl,
                                                  N=N)

        if min_cost is not None:
            if min_cost < 0:
                min_cost = 0
            profit_attack_with_m = EV_fork_markov_blocks_effort(p_B=p_B,
                                                          p_E=p_E,
                                                          p_V=p_V,
                                                          p_I=p_I,
                                                          p_m=p_m,
                                                          z=z,
                                                          n_fork=n_fork,
                                                          n_main=n_main,
                                                          epsilon=min_cost,
                                                          kr=kr,
                                                          kl=kl,
                                                          N=N)

        else:
            min_cost = float("NaN")
            profit_attack_with_m = 0

        prob_success_attack_without_m = Pr_success(kr=kr,
                   kl=kl,
                   N=N,
                   z=z,
                   p_B=p_B,
                   p_A=p_A,
                   p_I=p_I,
                   p_E=p_E,
                   p_V=p_V + p_m)[Prslt]

        prob_success_attack_with_m = Pr_success(kr=kr,
                   kl=kl,
                   N=N,
                   z=z,
                   p_B=p_B,
                   p_A=p_A,
                   p_I=p_I,
                   p_E=p_E + p_m,
                   p_V=p_V)[Prslt]
        """
        if profit_attack_with_m > proift_main_with_m:
            mark = "X"
        else:
            mark = " "
        print("p_E = {:4.2f}  p_A = {:4.2f}  p_m = {:4.2f}  EV_main = {:.3f}  P_abstain = {:.2f}  e = {:.2f}  EV_fork = {:.3f}  P_join = {:.3f} {}".format(
                                        p_E,
                                        1 - (p_E + p_m),
                                        p_m,
                                        proift_main_with_m,
                                        1 - prob_success_attack_without_m,
                                        min_cost,
                                        profit_attack_with_m,
                                        prob_success_attack_with_m,
                                        mark))
        """
        if min_cost <= 0.0:
            table += "\t& \\cellcolor{blue!15}"
        else:
            table += "\t& "
        table += """\\makecell[l]{{ $p_\\mathcal{{A}} = {:.3f} $ \\\\ $\\rho= {:.3f}$ \\\\ $\\epsilon = {:.3f}$ \\\\ $\\rho'= {:.3f}$ \\\\ $\\Pr = {:.3f}$ }} """.format(
                        1 - (p_B + p_E + p_V + p_m + p_I),
                        proift_main_with_m,
                        min_cost,
                        profit_attack_with_m,
                        prob_success_attack_with_m)

table +="""
\\end{{tabular}}
}}
\\caption{{Comparison of minimum bribe required per block $ \\epsilon $ for $\\rho_{{\\textit{{fork comp.}}}}>\\rho_{{\\textit{{main comp.}}}}$ with effort-related compensation and probabilities calculated using our Markov chain.
The axis iterate the hashrate of an individual miner $ p_m $, and other attackers $ p_\\mathcal{{E}} $.
The table also shows the expected reward of miner $ m $, if $ p_m $ would be directed towards the attack chain $                \\rho'=\\rho_{{\\textit{{fork comp.}}}} $, as well as the expected reward $ \\rho=\\rho_{{\\textit{{main comp.}}}} $, if $ p_m $    would be directed towards the main chain.
All attacks start with a disadvantage of $ \\kbe = {:d} $ and a duration of $ N = {:d} $ and the following configuration of the Markov chain $\\dkra = {:d}, \\dkla = {:d}, \\eta_{{attack}} = {:d}, \\eta_{{main}} = {:d}, p_\\mathcal{{V}} = {:f}$.}}
\\label{{tab:txexclib}}
\\end{{table*}}
""".format(z,N,kr,kl,n_fork,n_main,p_V)

print(table)
#%store table >../paper/preprint/texsrc/tbl_profits_exclusion_markov_0.tex
