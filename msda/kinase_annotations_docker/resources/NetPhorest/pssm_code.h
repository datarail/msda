o = pssm(s, scansite_human_BRCT_BRCA1_pssm, 15)/2.258e+05;
if (o > 0) {
  o = log(o);
  o = 2.31706127099737e-07+(0.210947016967266-2.31706127099737e-07)/(1+exp(3.04234*(-5.94869-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tBRCT\tany_group\t%.6f\t%.6f\n", o, 0.0699854212223765);
}

o = pssm(s, scansite_human_KIN_AMPK_group_pssm, 15)/4.536e+04;
if (o > 0) {
  o = log(o);
  o = 0.00188008701167695+(0.106172644182624-0.00188008701167695)/(1+exp(0.598562*(-3.67134-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tAMPK_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
}

o = pssm(s, scansite_human_KIN_CDK5_pssm, 15)/1.034e+05;
if (o > 0) {
  o = log(o);
  o = 2.75740855529468e-08+(0.360744375400965-2.75740855529468e-08)/(1+exp(0.97807*(-3.18552-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tCDK2_CDK3_CDK1_CDK5_group\t%.6f\t%.6f\n", o, 0.04);
}

o = pssm(s, scansite_human_KIN_CK1_group_pssm, 15)/1.399e+05;
if (o > 0) {
  o = log(o);
  o = 9.95701610076016e-07+(0.145180988444129-9.95701610076016e-07)/(1+exp(0.993866*(-4.0268-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tCK1_group\t%.6f\t%.6f\n", o, 0.0529150262212918);
}

o = pssm(s, scansite_human_KIN_DMPK1_pssm, 15)/2.688e+03;
if (o > 0) {
  o = log(o);
  o = 0.00918460458511681+(0.291871921182266-0.00918460458511681)/(1+exp(0.63519*(-1.97074-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tDMPK_group\t%.6f\t%.6f\n", o, 0.04);
}

o = pssm(s, scansite_human_KIN_EGFR_pssm, 15)/9.997e+05;
if (o > 0) {
  o = log(o);
  o = 0.0196920514889329+(0.0813793809332912-0.0196920514889329)/(1+exp(3.73875*(-3.21577-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tEGFR_group\t%.6f\t%.6f\n", o, 0.04);
}

o = pssm(s, scansite_human_KIN_InsR_pssm, 15)/2.190e+04;
if (o > 0) {
  o = log(o);
  o = 0.0237693672837333+(0.102735436709296-0.0237693672837333)/(1+exp(3.09673*(-2.53617-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tInsR_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
}

o = pssm(s, scansite_human_KIN_Lck_pssm, 15)/1.877e+05;
if (o > 0) {
  o = log(o);
  o = 0.0237856170384623+(0.151954202989566-0.0237856170384623)/(1+exp(0.999999*(-3-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tSrc_group\t%.6f\t%.6f\n", o, 0.066332495807108);
}

o = pssm(s, scansite_human_KIN_MAPK14_pssm, 15)/1.485e+04;
if (o > 0) {
  o = log(o);
  o = 1.19929673242453e-08+(0.148220256486223-1.19929673242453e-08)/(1+exp(43.6605*(-2.82847-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tp38_group\t%.6f\t%.6f\n", o, 0.04);
}

o = pssm(s, scansite_human_KIN_MAPK3_pssm, 15)/2.383e+05;
if (o > 0) {
  o = log(o);
  o = 0.000500417034004853+(0.331151825906694-0.000500417034004853)/(1+exp(0.659862*(-2.14883-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tMAPK3_MAPK1_MAPK7_NLK_group\t%.6f\t%.6f\n", o, 0.04);
}

o = pssm(s, scansite_human_KIN_PKD1_pssm, 15)/1.814e+06;
if (o > 0) {
  o = log(o);
  o = 3.91054175771969e-08+(0.116863312249919-3.91054175771969e-08)/(1+exp(1.64269*(-4.77718-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tKIN\tPKD_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
}

o = pssm(s, scansite_human_PTB_SHC1_pssm, 15)/2.583e+04;
if (o > 0) {
  o = log(o);
  o = 1.44098589574965e-08+(0.577823050045167-1.44098589574965e-08)/(1+exp(1.87041*(-1.73069-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tPTB\tSHC1_SHC2_SHC3_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
}

o = pssm(s, scansite_human_SH2_PLCG1_1_pssm, 15)/4.975e+04;
if (o > 0) {
  o = log(o);
  o = 0.00804220390320895+(0.131040568715424-0.00804220390320895)/(1+exp(1.1319*(-2.14053-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tSH2\tPLCG_1_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
}

o = pssm(s, scansite_human_SH2_PLCG1_2_pssm, 15)/5.961e+04;
if (o > 0) {
  o = log(o);
  o = 0.0285162648018724+(0.179889183290924-0.0285162648018724)/(1+exp(0.999999*(-2-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite\thuman\tSH2\tPLCG_2_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_INPP5D_pssm, 15)/1.387e-03;
if (o > 0) {
  o = log(o);
  o = 0.0122117256875409+(0.0823386607323488-0.0122117256875409)/(1+exp(1.27666*(-3.09424-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tINPP5D\t%.6f\t%.6f\n", o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_SH3BP2_pssm, 15)/2.135e-05;
if (o > 0) {
  o = log(o);
  o = 0.017911233061943+(0.116550893235419-0.017911233061943)/(1+exp(3.84889*(-2.44766-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSH3BP2\t%.6f\t%.6f\n", o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_SHC2_pssm, 15)/1.528e-01;
if (o > 0) {
  o = log(o);
  o = 0.0207364427656338+(0.101906853305672-0.0207364427656338)/(1+exp(1.31642*(-2.20547-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSHC2_SHC3_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_SHE_pssm, 15)/2.430e-03;
if (o > 0) {
  o = log(o);
  o = 0.0116298990401533+(0.133810106332926-0.0116298990401533)/(1+exp(0.608993*(-1.51017-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSHE\t%.6f\t%.6f\n", o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_Syk_1_pssm, 15)/4.793e-03;
if (o > 0) {
  o = log(o);
  o = 0.00656765550958794+(0.10175596578118-0.00656765550958794)/(1+exp(1.527*(-3.1099-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tSyk_1\t%.6f\t%.6f\n", o, 0.0285714285714286);
}

o = pssm(s, scansite_LiLabOPAL_human_SH2_TNS4_pssm, 15)/1.657e-03;
if (o > 0) {
  o = log(o);
  o = 0.00236464040993422+(0.137621910237852-0.00236464040993422)/(1+exp(0.75074*(-2.58628-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_LiLabOPAL\thuman\tSH2\tTNS4\t%.6f\t%.6f\n", o, 0.0285714285714286);
}


o = pssm(s, scansite_TurkLabOPL_human_KIN_CLK2_pssm, 15)/1.222e+03;
if (o > 0) {
  o = log(o);
  o = 0.00503259201635941+(0.421209858103062-0.00503259201635941)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tCLK_group\t%.6f\t%.6f\n", o, 0.04);
}


o = pssm(s, scansite_TurkLabOPL_human_KIN_DAPK3_pssm, 15)/6.181e+02;
if (o > 0) {
  o = log(o);
  o = 9.25307706661765e-11+(0.112107567281039-9.25307706661765e-11)/(1+exp(99.4089*(-13.4476-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tDAPK_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
}


o = pssm(s, scansite_TurkLabOPL_human_KIN_ICK_pssm, 15)/1.632e+03;
if (o > 0) {
  o = log(o);
  o = 0.00436129308904892+(0.386594515151589-0.00436129308904892)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tRCK_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_LKB1_pssm, 15)/3.898e+02;
if (o > 0) {
  o = log(o);
  o = 0.000122769270670808+(0.110169491525424-0.000122769270670808)/(1+exp(4.61123*(-2.49789-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tLKB1\t%.6f\t%.6f\n", o, 0.02);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_MST1_pssm, 15)/1.817e+03;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tMST_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_MST4_pssm, 15)/1.521e+03;
if (o > 0) {
  o = log(o);
  o = 0.00436129308904892+(0.386594515151589-0.00436129308904892)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tYSK_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_PAK2_pssm, 15)/6.245e+02;
if (o > 0) {
  o = log(o);
  o = 0.010462108110282+(0.190793697897193-0.010462108110282)/(1+exp(1.95923*(-3.42444-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tPAK_group\t%.6f\t%.6f\n", o, 0.0489897948556636);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_Pim1_pssm, 15)/2.386e+03;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tPim3_Pim1_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_Pim2_pssm, 15)/1.843e+03;
if (o > 0) {
  o = log(o);
  o = 0.0025226437265188+(0.266792809839167-0.0025226437265188)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tPim2\t%.6f\t%.6f\n", o, 0.02);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_SLK_pssm, 15)/5.006e+02;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tSLK_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_TGFbR2_pssm, 15)/6.444e+02;
if (o > 0) {
  o = log(o);
  o = 0.00436129308904892+(0.386594515151589-0.00436129308904892)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tACTR2_ACTR2B_TGFbR2_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_TLK1_pssm, 15)/1.819e+03;
if (o > 0) {
  o = log(o);
  o = 0.00356383307469715+(0.339755861506917-0.00356383307469715)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tTLK_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL_human_KIN_TNIK_pssm, 15)/1.566e+03;
if (o > 0) {
  o = log(o);
  o = 0.00503259201635941+(0.421209858103062-0.00503259201635941)/(1+exp(0.617646*(-1.89176-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tMSN_group\t%.6f\t%.6f\n", o, 0.04);
}


o = pssm(s, scansite_TurkLabOPL_human_KIN_p70S6K_pssm, 15)/1.776e+03;
if (o > 0) {
  o = log(o);
  o = 0.00412223670003529+(0.0660059353402013-0.00412223670003529)/(1+exp(1.06901*(-3.53051-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL\thuman\tKIN\tp70S6K_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
}

o = pssm(s, scansite_TurkLabOPL1_human_KIN_EphA3_pssm, 15)/3.156e+02;
if (o > 0) {
  o = log(o);
  o = 0.048195566335724+(0.144279729622026-0.048195566335724)/(1+exp(8.36245*(-1.70993-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_TurkLabOPL1\thuman\tKIN\tEph_group\t%.6f\t%.6f\n", o, 0.0748331477354788);
}



o = pssm(s, scansite_YaffeLabOPL_human_KIN_NEK2_pssm, 15)/9.353e+02;
if (o > 0) {
  o = log(o);
  o = 0.0296541137013766+(0.144003522831755-0.0296541137013766)/(1+exp(2.00528*(-1.93177-o)));
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tscansite_YaffeLabOPL\thuman\tKIN\tNEK1_NEK5_NEK3_NEK4_NEK11_NEK2_group\t%.6f\t%.6f\n", o, 0.0489897948556636);
}



