"""Provisional research figure; exact data and correlated walk estimates separated."""
try:  # put every scripts/ topic folder on the import path
    import _scripts_path  # noqa: F401
except ImportError:  # running from a flat copy of the scripts
    pass
import argparse
import json
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser();p.add_argument('results',type=Path);args=p.parse_args()
    r=json.loads(args.results.read_text());out=args.results.parent
    with plt.rc_context({'font.size':10,'pdf.fonttype':42,'figure.facecolor':'white'}):
        fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
        ax=axes[0,0]
        for low,label,color,mark in ((0,'<1,000 states','#333333','o'),
                                     (1000,'1,000–9,999 states','#0072B2','s'),
                                     (10000,'≥10,000 states','#A34E00','^')):
            z=[b for b in r['exact']['size_bins_by_genus'] if b['size_low']==low]
            ax.plot([b['genus'] for b in z],[b['median_tv'] for b in z],marker=mark,color=color,label=label)
        ax.set(xlabel='Genus',ylabel='Median component TV distance',ylim=(0,.7),xticks=range(1,6),
               title='A. Exact S₄–S₇ components: size matters')
        ax.legend(frameon=False,fontsize=8)
        ax=axes[0,1]
        for a,label,color,marker in (('S','Uniform S_d starts','#0072B2','o'),('A','Uniform A_d starts','#A34E00','s')):
            z=[max((x for x in r['iid'] if x['degree']==d and x['ambient']==a),key=lambda x:x['n']) for d in (8,9,12,18,25)]
            ax.plot([x['degree'] for x in z],[x['tv_coin'] for x in z],marker=marker,color=color,label=label)
        ax.set(xlabel='Degree d',ylabel='TV distance from genus-coin model',ylim=(0,.55),xticks=[8,12,18,25],
               title='B. Independent starts: most frequent genus')
        ax.legend(frameon=False,fontsize=8)
        ax.text(.03,.96,'Genus differs between degrees; ≥76,000 retained starts/cell',transform=ax.transAxes,va='top',fontsize=8)
        ax=axes[1,0]
        # No walk error bars: observations are dependent, mixing unproved.
        for g in range(1,8):
            z=[w for w in r['walks'] if w['genus']==g]
            ax.scatter([g]*len(z),[w['tv_coin'] for w in z],s=35,facecolors='none',edgecolors='#333333')
            ax.plot(g,median(w['tv_coin'] for w in z),'_',color='#0072B2',markersize=18)
        ax.set(xlabel='Genus',ylabel='Estimated component TV distance',ylim=(0,.9),xticks=range(1,8),
               title='C. Selected S₈ components: walk estimates')
        ax.text(.02,.98,'5 selected components/genus; blue = median\n128,000 correlated observations/component',transform=ax.transAxes,va='top',fontsize=8)
        ax=axes[1,1]
        x=next(x for x in r['iid'] if x['degree']==25 and x['ambient']=='A' and x['genus']==19)
        import sys
        sys.path.insert(0,str(Path(__file__).resolve().parent))
        from defect_coin_followup import coin_float
        k=list(range(9));degree=coin_float(24)
        degree=[2*p if j%2==0 else 0 for j,p in enumerate(degree)]
        ax.plot(k,x['observed'][:9],color='#333333',marker='o',label='Observed A₂₅ starts')
        ax.plot(k,x['coin'][:9],color='#A34E00',marker='s',linestyle='--',label='Genus 19 coins')
        ax.plot(k,degree[:9],color='#0072B2',marker='^',linestyle=':',label='Degree cycles, even parity')
        ax.set(xlabel='Defect δ',ylabel='Probability',ylim=(0,.55),xticks=range(9),
               title='D. Parity explains the alternating mismatch')
        ax.legend(frameon=False,fontsize=8)
        for ax in axes.flat:
            ax.spines[['top','right']].set_visible(False)
            ax.grid(axis='y',alpha=.2)
        fig.suptitle('Defect versus uneven coins: exact census and exploratory sampling',fontsize=13)
        for ext in ('png','pdf'):
            fig.savefig(out/f'COMPARISON.{ext}',dpi=170,facecolor='white')
    print(out/'COMPARISON.png')


if __name__=='__main__':main()
