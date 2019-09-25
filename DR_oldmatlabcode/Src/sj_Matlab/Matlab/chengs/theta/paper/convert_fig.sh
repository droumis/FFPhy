dvips -i -S1 figs -o Cheng_fig
Smv "" ps Cheng_fig00?

mv Cheng_fig001.ps v3.3/Cheng_Fig1.ps
mv Cheng_fig002.ps v3.3/Cheng_Fig2.ps
mv Cheng_fig003.ps v3.3/Cheng_Fig3.ps
mv Cheng_fig004.ps v3.3/Cheng_Fig4.ps

Sps2eps Cheng_fig001
Sps2eps Cheng_fig002
Sps2eps Cheng_fig003
Sps2eps Cheng_fig004

convert Cheng_fig001.eps Cheng_fig001.jpg
convert Cheng_fig002.eps Cheng_fig002.jpg
convert Cheng_fig003.eps Cheng_fig003.jpg
convert Cheng_fig004.eps Cheng_fig004.jpg
