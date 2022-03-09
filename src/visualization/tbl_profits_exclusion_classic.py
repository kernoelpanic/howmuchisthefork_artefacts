#!/usr/bin/env python3

from profitability import *

z=1
# init all variables to zero
n_main=0
n_fork=0
p_B=0
p_A=None
p_I=0
p_V=0
p_E=0
p_m=0

table = """
\\renewcommand{\\arraystretch}{1.3} % increase space between lines
\\begin{table*}[!htpb]
\\centering
\\scriptsize
\\resizebox{1.0\\linewidth}{!}{   % 0.7 for 2 colum view
\\begin{tabular}{c|c|c|c|c|c|c|c|c|}
"""

O_m = [0.05, 0.10, 0.20, 0.30, 0.33, 0.382, 0.40, ]
O_E = [0.0, 0.05, 0.10, 0.20, 0.30, 0.33, 0.382, 0.40, ]

for i in O_m:
    table += "\t& $ p_m = $ {} ".format(str(i))
#print("\\ \\hline")

for p_E in O_E:
    table += """\\\\ \\hline
    """
    table += "$ p_\\mathcal{{ E }} = {:.2f} $".format(p_E)
    for p_m in O_m:
        if 1 - (p_E + p_m) < 0:
            continue

        proift_main_with_m = EV_main_blocks(p_B=p_B,
                                           p_E=p_E,
                                           p_V=p_V,
                                           p_I=p_I,
                                           p_m=p_m,
                                           z=z,
                                           n_main=n_main)

        min_cost = e_min(p_B=p_B,
                                        p_E=p_E,
                                        p_V=p_V,
                                        p_I=p_I,
                                        p_m=p_m,
                                        z=z)

        if min_cost is not None:
            if min_cost < 0:
                min_cost = 0
            profit_attack_with_m = EV_fork_blocks(p_B=p_B,
                                                   p_E=p_E,
                                                   p_V=p_V,
                                                   p_I=p_I,
                                                   p_m=p_m,
                                                   z=z,
                                                   n_fork=n_fork,
                                                   epsilon=min_cost)
        else:
            min_cost = float("NaN")
            profit_attack_with_m = 0

        prob_success_attack_without_m = Pr_success_abstain(p_E=p_E,
                                                           p_m=p_m)
        prob_success_attack_with_m = Pr_success_join(p_E=p_E,
                                                     p_m=p_m)
        """
        if profit_attack_with_m > proift_main_with_m:
            mark = "X"
        else:
            mark = " "
        print("p_E = {:4.2f}  p_A = {:4.2f}  p_m = {:4.2f}  EV_main = {:.3f}  Pr_abstain = {:.2f}  e = {:.2f}  EV_fork = {:.3f}  Pr_join = {:.3f} {}".format(
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

table += """
\\end{tabular}
}
\\caption{Comparison of minimum bribe required per block $ \\epsilon $ for $\\rho_{fork}>\\rho_{main}$ with classic infinite probability calculation.
The axis iterate the hashrate of an individual miner $ p_m $, and other attackers $ p_\\mathcal{E} $.
The table also shows the expected reward of miner $ m $, if $ p_m $ would be directed towards the attack chain $ \\rho'=\\rho_{fork} $, as well as the expected reward $ \\rho=\\rho_{main} $, if $ p_m $ would be directed towards the main chain.
All attacks start with a disadvantage of $ \\kbe = 1 $.}
\\label{tab:txexclib}
\\end{table*}
"""

print(table)

# Notebook extract string to file:
#%store table >../paper/preprint/texsrc/tbl_profits_exclusion.tex
