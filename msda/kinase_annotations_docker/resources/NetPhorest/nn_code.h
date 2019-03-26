if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_1, 7, 15);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_2, 5, 20);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_3, 7, 2);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_4, 5, 20);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_5, 13, 10);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_6, 11, 2);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_7, 5, 20);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_8, 7, 20);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_9, 9, 0);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_10, 5, 4);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_11, 7, 15);
  o += feed_forward(s, nn_human_KIN_ATM_ATR_group_12, 9, 0);
  o /= 12;
  o = 8.14564048016133e-07+(0.464280740430288-8.14564048016133e-07)/(1+exp(84.3464*(0.496148-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tATM_ATR_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_Abl_group_1, 5, 10);
  o += feed_forward(s, nn_human_KIN_Abl_group_2, 9, 0);
  o += feed_forward(s, nn_human_KIN_Abl_group_3, 7, 6);
  o += feed_forward(s, nn_human_KIN_Abl_group_4, 7, 4);
  o += feed_forward(s, nn_human_KIN_Abl_group_5, 5, 10);
  o += feed_forward(s, nn_human_KIN_Abl_group_6, 11, 6);
  o += feed_forward(s, nn_human_KIN_Abl_group_7, 7, 15);
  o += feed_forward(s, nn_human_KIN_Abl_group_8, 11, 6);
  o += feed_forward(s, nn_human_KIN_Abl_group_9, 7, 4);
  o += feed_forward(s, nn_human_KIN_Abl_group_10, 9, 4);
  o += feed_forward(s, nn_human_KIN_Abl_group_11, 7, 4);
  o += feed_forward(s, nn_human_KIN_Abl_group_12, 7, 20);
  o /= 12;
  o = 0.0220853111848124+(0.0754387553643689-0.0220853111848124)/(1+exp(99.9989*(0.503018-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tAbl_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_AuroraA_1, 13, 20);
  o += feed_forward(s, nn_human_KIN_AuroraA_2, 9, 0);
  o += feed_forward(s, nn_human_KIN_AuroraA_3, 9, 2);
  o += feed_forward(s, nn_human_KIN_AuroraA_4, 11, 10);
  o += feed_forward(s, nn_human_KIN_AuroraA_5, 5, 6);
  o += feed_forward(s, nn_human_KIN_AuroraA_6, 7, 2);
  o += feed_forward(s, nn_human_KIN_AuroraA_7, 11, 20);
  o += feed_forward(s, nn_human_KIN_AuroraA_8, 5, 10);
  o += feed_forward(s, nn_human_KIN_AuroraA_9, 5, 10);
  o += feed_forward(s, nn_human_KIN_AuroraA_10, 9, 20);
  o += feed_forward(s, nn_human_KIN_AuroraA_11, 7, 15);
  o += feed_forward(s, nn_human_KIN_AuroraA_12, 7, 2);
  o /= 12;
  o = 0.00827755294440268+(0.0463703229932843-0.00827755294440268)/(1+exp(99.9999*(0.485745-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tAuroraA\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_1, 7, 2);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_2, 5, 0);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_3, 7, 0);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_4, 13, 0);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_5, 5, 2);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_6, 5, 0);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_7, 7, 2);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_8, 5, 0);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_9, 5, 6);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_10, 13, 10);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_11, 7, 2);
  o += feed_forward(s, nn_human_KIN_CDK4_CDK6_group_12, 5, 0);
  o /= 12;
  o = 0.0025408733878434+(0.125418782750793-0.0025408733878434)/(1+exp(67.4791*(0.489029-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tCDK4_CDK6_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_CDK7_1, 5, 20);
  o += feed_forward(s, nn_human_KIN_CDK7_2, 13, 20);
  o += feed_forward(s, nn_human_KIN_CDK7_3, 9, 2);
  o += feed_forward(s, nn_human_KIN_CDK7_4, 5, 20);
  o += feed_forward(s, nn_human_KIN_CDK7_5, 13, 20);
  o += feed_forward(s, nn_human_KIN_CDK7_6, 5, 6);
  o += feed_forward(s, nn_human_KIN_CDK7_7, 5, 2);
  o += feed_forward(s, nn_human_KIN_CDK7_8, 7, 4);
  o += feed_forward(s, nn_human_KIN_CDK7_9, 9, 2);
  o += feed_forward(s, nn_human_KIN_CDK7_10, 5, 0);
  o += feed_forward(s, nn_human_KIN_CDK7_11, 13, 20);
  o += feed_forward(s, nn_human_KIN_CDK7_12, 13, 20);
  o /= 12;
  o = 0.00880910147557574+(0.0321862948679917-0.00880910147557574)/(1+exp(65.1908*(0.504986-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tCDK7\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_CK2_group_1, 11, 15);
  o += feed_forward(s, nn_human_KIN_CK2_group_2, 7, 15);
  o += feed_forward(s, nn_human_KIN_CK2_group_3, 9, 15);
  o += feed_forward(s, nn_human_KIN_CK2_group_4, 9, 20);
  o += feed_forward(s, nn_human_KIN_CK2_group_5, 11, 20);
  o += feed_forward(s, nn_human_KIN_CK2_group_6, 7, 20);
  o += feed_forward(s, nn_human_KIN_CK2_group_7, 11, 10);
  o += feed_forward(s, nn_human_KIN_CK2_group_8, 9, 20);
  o += feed_forward(s, nn_human_KIN_CK2_group_9, 13, 0);
  o += feed_forward(s, nn_human_KIN_CK2_group_10, 9, 0);
  o += feed_forward(s, nn_human_KIN_CK2_group_11, 13, 2);
  o += feed_forward(s, nn_human_KIN_CK2_group_12, 9, 10);
  o /= 12;
  o = 0.00108607534568641+(0.556057610497607-0.00108607534568641)/(1+exp(19.3479*(0.543655-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tCK2_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_CaMKII_group_1, 7, 0);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_2, 7, 6);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_3, 13, 15);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_4, 13, 2);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_5, 9, 2);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_6, 13, 10);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_7, 13, 15);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_8, 11, 20);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_9, 7, 0);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_10, 7, 15);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_11, 7, 6);
  o += feed_forward(s, nn_human_KIN_CaMKII_group_12, 13, 10);
  o /= 12;
  o = 0.0109588617056442+(0.165880291183607-0.0109588617056442)/(1+exp(36.9773*(0.494663-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tCaMKII_group\t%.6f\t%.6f\n", o, 0.04);
  }
}

if (c == 'S') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_CaMKIV_1, 5, 10);
  o += feed_forward(s, nn_human_KIN_CaMKIV_2, 11, 4);
  o += feed_forward(s, nn_human_KIN_CaMKIV_3, 7, 10);
  o += feed_forward(s, nn_human_KIN_CaMKIV_4, 5, 6);
  o += feed_forward(s, nn_human_KIN_CaMKIV_5, 11, 4);
  o += feed_forward(s, nn_human_KIN_CaMKIV_6, 11, 4);
  o += feed_forward(s, nn_human_KIN_CaMKIV_7, 7, 10);
  o += feed_forward(s, nn_human_KIN_CaMKIV_8, 7, 20);
  o += feed_forward(s, nn_human_KIN_CaMKIV_9, 9, 20);
  o /= 9;
  o = 0.00706832716257153+(0.0441643719550943-0.00706832716257153)/(1+exp(100*(0.499798-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tCaMKIV\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_CaMKI_group_1, 13, 4);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_2, 5, 20);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_3, 5, 4);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_4, 7, 15);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_5, 5, 20);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_6, 7, 20);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_7, 7, 15);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_8, 7, 20);
  o += feed_forward(s, nn_human_KIN_CaMKI_group_9, 5, 20);
  o /= 9;
  o = 0.0148247940170952+(0.068663257852447-0.0148247940170952)/(1+exp(99.9999*(0.508282-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tCaMKI_group\t%.6f\t%.6f\n", o, 0.04);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_DNAPK_1, 7, 0);
  o += feed_forward(s, nn_human_KIN_DNAPK_2, 5, 4);
  o += feed_forward(s, nn_human_KIN_DNAPK_3, 7, 15);
  o += feed_forward(s, nn_human_KIN_DNAPK_4, 11, 15);
  o += feed_forward(s, nn_human_KIN_DNAPK_5, 9, 0);
  o += feed_forward(s, nn_human_KIN_DNAPK_6, 7, 15);
  o += feed_forward(s, nn_human_KIN_DNAPK_7, 5, 0);
  o += feed_forward(s, nn_human_KIN_DNAPK_8, 7, 15);
  o += feed_forward(s, nn_human_KIN_DNAPK_9, 7, 15);
  o += feed_forward(s, nn_human_KIN_DNAPK_10, 5, 0);
  o += feed_forward(s, nn_human_KIN_DNAPK_11, 5, 15);
  o += feed_forward(s, nn_human_KIN_DNAPK_12, 7, 0);
  o /= 12;
  o = 3.96139289555181e-07+(0.129453981998291-3.96139289555181e-07)/(1+exp(58.6663*(0.476378-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tDNAPK\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_EIF2AK2_1, 5, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_2, 9, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_3, 9, 15);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_4, 9, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_5, 9, 15);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_6, 5, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_7, 5, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_8, 9, 15);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_9, 5, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_10, 9, 2);
  o += feed_forward(s, nn_human_KIN_EIF2AK2_11, 5, 0);
  o /= 11;
  o = 0.00864056979850585+(0.0431188796462041-0.00864056979850585)/(1+exp(99.9999*(0.499404-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tEIF2AK2\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_1, 9, 6);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_2, 5, 2);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_3, 13, 15);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_4, 7, 2);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_5, 9, 6);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_6, 5, 2);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_7, 5, 2);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_8, 7, 15);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_9, 5, 10);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_10, 7, 10);
  o += feed_forward(s, nn_human_KIN_FLT3_CSF1R_Kit_PDGFR_group_11, 5, 2);
  o /= 11;
  o = 0.0128617281449406+(0.10844004778682-0.0128617281449406)/(1+exp(100*(0.495483-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tFLT3_CSF1R_Kit_PDGFR_group\t%.6f\t%.6f\n", o, 0.0447213595499958);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_GRK_group_1, 9, 4);
  o += feed_forward(s, nn_human_KIN_GRK_group_2, 9, 10);
  o += feed_forward(s, nn_human_KIN_GRK_group_3, 9, 6);
  o += feed_forward(s, nn_human_KIN_GRK_group_4, 7, 20);
  o += feed_forward(s, nn_human_KIN_GRK_group_5, 5, 10);
  o += feed_forward(s, nn_human_KIN_GRK_group_6, 11, 15);
  o += feed_forward(s, nn_human_KIN_GRK_group_7, 13, 0);
  o += feed_forward(s, nn_human_KIN_GRK_group_8, 13, 4);
  o += feed_forward(s, nn_human_KIN_GRK_group_9, 7, 20);
  o += feed_forward(s, nn_human_KIN_GRK_group_10, 9, 20);
  o += feed_forward(s, nn_human_KIN_GRK_group_11, 7, 6);
  o += feed_forward(s, nn_human_KIN_GRK_group_12, 13, 4);
  o /= 12;
  o = 0.0128877078967273+(0.1168316364453-0.0128877078967273)/(1+exp(99.8727*(0.417367-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tGRK_group\t%.6f\t%.6f\n", o, 0.0529150262212918);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_GSK3_group_1, 13, 4);
  o += feed_forward(s, nn_human_KIN_GSK3_group_2, 9, 15);
  o += feed_forward(s, nn_human_KIN_GSK3_group_3, 9, 20);
  o += feed_forward(s, nn_human_KIN_GSK3_group_4, 13, 0);
  o += feed_forward(s, nn_human_KIN_GSK3_group_5, 9, 15);
  o += feed_forward(s, nn_human_KIN_GSK3_group_6, 9, 4);
  o += feed_forward(s, nn_human_KIN_GSK3_group_7, 9, 15);
  o += feed_forward(s, nn_human_KIN_GSK3_group_8, 9, 4);
  o += feed_forward(s, nn_human_KIN_GSK3_group_9, 13, 4);
  o += feed_forward(s, nn_human_KIN_GSK3_group_10, 13, 15);
  o += feed_forward(s, nn_human_KIN_GSK3_group_11, 11, 10);
  o += feed_forward(s, nn_human_KIN_GSK3_group_12, 13, 4);
  o /= 12;
  o = 1.48346095777565e-06+(0.15039430853318-1.48346095777565e-06)/(1+exp(15.9569*(0.470278-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tGSK3_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_1, 7, 0);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_2, 13, 0);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_3, 13, 0);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_4, 5, 4);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_5, 11, 6);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_6, 5, 0);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_7, 7, 4);
  o += feed_forward(s, nn_human_KIN_HIPK1_HIPK2_group_8, 5, 0);
  o /= 8;
  o = 0.000917679784190824+(0.0781078627402849-0.000917679784190824)/(1+exp(94.4907*(0.480654-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tHIPK1_HIPK2_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_1, 5, 4);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_2, 9, 15);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_3, 7, 15);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_4, 7, 10);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_5, 9, 2);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_6, 5, 4);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_7, 11, 10);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_8, 7, 20);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_9, 7, 15);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_10, 13, 10);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_11, 11, 10);
  o += feed_forward(s, nn_human_KIN_IKKalpha_IKKbeta_group_12, 7, 15);
  o /= 12;
  o = 0.00793378458421489+(0.059909298866143-0.00793378458421489)/(1+exp(99.9997*(0.498356-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tIKKalpha_IKKbeta_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_JAK2_1, 11, 10);
  o += feed_forward(s, nn_human_KIN_JAK2_2, 5, 4);
  o += feed_forward(s, nn_human_KIN_JAK2_3, 11, 4);
  o += feed_forward(s, nn_human_KIN_JAK2_4, 7, 10);
  o += feed_forward(s, nn_human_KIN_JAK2_5, 5, 0);
  o += feed_forward(s, nn_human_KIN_JAK2_6, 11, 4);
  o += feed_forward(s, nn_human_KIN_JAK2_7, 11, 0);
  o += feed_forward(s, nn_human_KIN_JAK2_8, 13, 20);
  o += feed_forward(s, nn_human_KIN_JAK2_9, 11, 4);
  o += feed_forward(s, nn_human_KIN_JAK2_10, 7, 4);
  o += feed_forward(s, nn_human_KIN_JAK2_11, 11, 4);
  o += feed_forward(s, nn_human_KIN_JAK2_12, 7, 2);
  o /= 12;
  o = 1.26323212117766e-08+(0.0514651588561949-1.26323212117766e-08)/(1+exp(71.3349*(0.481616-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tJAK2\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_JNK_group_1, 7, 4);
  o += feed_forward(s, nn_human_KIN_JNK_group_2, 13, 10);
  o += feed_forward(s, nn_human_KIN_JNK_group_3, 5, 6);
  o += feed_forward(s, nn_human_KIN_JNK_group_4, 5, 0);
  o += feed_forward(s, nn_human_KIN_JNK_group_5, 5, 2);
  o += feed_forward(s, nn_human_KIN_JNK_group_6, 5, 10);
  o += feed_forward(s, nn_human_KIN_JNK_group_7, 5, 6);
  o += feed_forward(s, nn_human_KIN_JNK_group_8, 9, 10);
  o += feed_forward(s, nn_human_KIN_JNK_group_9, 13, 0);
  o += feed_forward(s, nn_human_KIN_JNK_group_10, 7, 10);
  o += feed_forward(s, nn_human_KIN_JNK_group_11, 7, 0);
  o += feed_forward(s, nn_human_KIN_JNK_group_12, 9, 2);
  o /= 12;
  o = 0.00173002391707325+(0.139173807601364-0.00173002391707325)/(1+exp(69.0814*(0.492372-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tJNK_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_1, 5, 2);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_2, 5, 15);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_3, 9, 4);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_4, 9, 6);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_5, 9, 2);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_6, 5, 15);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_7, 13, 15);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_8, 9, 6);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_9, 5, 15);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_10, 5, 20);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_11, 5, 2);
  o += feed_forward(s, nn_human_KIN_KDR_FLT1_group_12, 5, 6);
  o /= 12;
  o = 0.00769461285325876+(0.0581995229190606-0.00769461285325876)/(1+exp(99.9986*(0.494094-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tKDR_FLT1_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T' || c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_MAP2K_group_1, 7, 10);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_2, 11, 10);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_3, 9, 2);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_4, 5, 15);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_5, 7, 0);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_6, 7, 4);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_7, 5, 2);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_8, 7, 10);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_9, 5, 10);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_10, 7, 10);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_11, 7, 2);
  o += feed_forward(s, nn_human_KIN_MAP2K_group_12, 11, 15);
  o /= 12;
  o = 0.0120436004400658+(0.170068227084073-0.0120436004400658)/(1+exp(72.9472*(0.496429-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tMAP2K_group\t%.6f\t%.6f\n", o, 0.0529150262212918);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_Met_group_1, 5, 10);
  o += feed_forward(s, nn_human_KIN_Met_group_2, 5, 4);
  o += feed_forward(s, nn_human_KIN_Met_group_3, 7, 4);
  o += feed_forward(s, nn_human_KIN_Met_group_4, 5, 4);
  o += feed_forward(s, nn_human_KIN_Met_group_5, 5, 6);
  o += feed_forward(s, nn_human_KIN_Met_group_6, 7, 4);
  o += feed_forward(s, nn_human_KIN_Met_group_7, 13, 6);
  o += feed_forward(s, nn_human_KIN_Met_group_8, 13, 2);
  o += feed_forward(s, nn_human_KIN_Met_group_9, 5, 10);
  o /= 9;
  o = 0.009966942951186+(0.0613358082975061-0.009966942951186)/(1+exp(100*(0.495183-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tMet_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_PDHK_group_1, 9, 10);
  o += feed_forward(s, nn_human_KIN_PDHK_group_2, 9, 20);
  o += feed_forward(s, nn_human_KIN_PDHK_group_3, 7, 0);
  o += feed_forward(s, nn_human_KIN_PDHK_group_4, 9, 2);
  o += feed_forward(s, nn_human_KIN_PDHK_group_5, 5, 0);
  o += feed_forward(s, nn_human_KIN_PDHK_group_6, 11, 15);
  o += feed_forward(s, nn_human_KIN_PDHK_group_7, 11, 6);
  o += feed_forward(s, nn_human_KIN_PDHK_group_8, 5, 2);
  o += feed_forward(s, nn_human_KIN_PDHK_group_9, 5, 20);
  o += feed_forward(s, nn_human_KIN_PDHK_group_10, 13, 4);
  o += feed_forward(s, nn_human_KIN_PDHK_group_11, 7, 6);
  o += feed_forward(s, nn_human_KIN_PDHK_group_12, 5, 0);
  o /= 12;
  o = 0.0183831573460718+(0.20795696526085-0.0183831573460718)/(1+exp(42.3627*(0.501103-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tPDHK_group\t%.6f\t%.6f\n", o, 0.0447213595499958);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_PKA_group_1, 13, 2);
  o += feed_forward(s, nn_human_KIN_PKA_group_2, 13, 10);
  o += feed_forward(s, nn_human_KIN_PKA_group_3, 9, 15);
  o += feed_forward(s, nn_human_KIN_PKA_group_4, 13, 20);
  o += feed_forward(s, nn_human_KIN_PKA_group_5, 11, 0);
  o += feed_forward(s, nn_human_KIN_PKA_group_6, 9, 20);
  o += feed_forward(s, nn_human_KIN_PKA_group_7, 13, 0);
  o += feed_forward(s, nn_human_KIN_PKA_group_8, 13, 10);
  o += feed_forward(s, nn_human_KIN_PKA_group_9, 9, 20);
  o += feed_forward(s, nn_human_KIN_PKA_group_10, 13, 6);
  o += feed_forward(s, nn_human_KIN_PKA_group_11, 9, 20);
  o += feed_forward(s, nn_human_KIN_PKA_group_12, 13, 2);
  o /= 12;
  o = 0.00264510988089579+(0.36196197658063-0.00264510988089579)/(1+exp(14.0501*(0.562018-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tPKA_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_PKB_group_1, 11, 10);
  o += feed_forward(s, nn_human_KIN_PKB_group_2, 11, 10);
  o += feed_forward(s, nn_human_KIN_PKB_group_3, 11, 15);
  o += feed_forward(s, nn_human_KIN_PKB_group_4, 11, 15);
  o += feed_forward(s, nn_human_KIN_PKB_group_5, 13, 10);
  o += feed_forward(s, nn_human_KIN_PKB_group_6, 11, 0);
  o += feed_forward(s, nn_human_KIN_PKB_group_7, 11, 2);
  o += feed_forward(s, nn_human_KIN_PKB_group_8, 11, 20);
  o += feed_forward(s, nn_human_KIN_PKB_group_9, 11, 20);
  o += feed_forward(s, nn_human_KIN_PKB_group_10, 11, 15);
  o += feed_forward(s, nn_human_KIN_PKB_group_11, 13, 10);
  o /= 11;
  o = 0.00196253321706196+(0.370717009236852-0.00196253321706196)/(1+exp(20.6761*(0.507834-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tPKB_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_PKC_group_1, 7, 6);
  o += feed_forward(s, nn_human_KIN_PKC_group_2, 13, 15);
  o += feed_forward(s, nn_human_KIN_PKC_group_3, 13, 2);
  o += feed_forward(s, nn_human_KIN_PKC_group_4, 5, 6);
  o += feed_forward(s, nn_human_KIN_PKC_group_5, 9, 6);
  o += feed_forward(s, nn_human_KIN_PKC_group_6, 7, 15);
  o += feed_forward(s, nn_human_KIN_PKC_group_7, 9, 4);
  o += feed_forward(s, nn_human_KIN_PKC_group_8, 11, 20);
  o += feed_forward(s, nn_human_KIN_PKC_group_9, 7, 15);
  o += feed_forward(s, nn_human_KIN_PKC_group_10, 13, 15);
  o += feed_forward(s, nn_human_KIN_PKC_group_11, 7, 20);
  o += feed_forward(s, nn_human_KIN_PKC_group_12, 9, 20);
  o /= 12;
  o = 3.27730558560227e-05+(0.525486467785397-3.27730558560227e-05)/(1+exp(30.1068*(0.520339-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tPKC_group\t%.6f\t%.6f\n", o, 0.06);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_1, 13, 4);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_2, 5, 6);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_3, 5, 20);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_4, 5, 2);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_5, 7, 10);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_6, 7, 2);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_7, 7, 15);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_8, 9, 0);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_9, 7, 10);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_10, 7, 10);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_11, 11, 15);
  o += feed_forward(s, nn_human_KIN_PKGcGK_group_12, 5, 6);
  o /= 12;
  o = 0.0198154247199193+(0.0702332672343646-0.0198154247199193)/(1+exp(100*(0.503187-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tPKGcGK_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_ROCK_group_1, 5, 2);
  o += feed_forward(s, nn_human_KIN_ROCK_group_2, 11, 20);
  o += feed_forward(s, nn_human_KIN_ROCK_group_3, 9, 10);
  o += feed_forward(s, nn_human_KIN_ROCK_group_4, 5, 15);
  o += feed_forward(s, nn_human_KIN_ROCK_group_5, 5, 4);
  o += feed_forward(s, nn_human_KIN_ROCK_group_6, 5, 2);
  o += feed_forward(s, nn_human_KIN_ROCK_group_7, 5, 2);
  o += feed_forward(s, nn_human_KIN_ROCK_group_8, 9, 15);
  o += feed_forward(s, nn_human_KIN_ROCK_group_9, 5, 20);
  o += feed_forward(s, nn_human_KIN_ROCK_group_10, 5, 10);
  o += feed_forward(s, nn_human_KIN_ROCK_group_11, 5, 2);
  o += feed_forward(s, nn_human_KIN_ROCK_group_12, 5, 15);
  o /= 12;
  o = 0.0179791879440675+(0.0774757891048534-0.0179791879440675)/(1+exp(99.9992*(0.513959-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tROCK_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_RSK_group_1, 11, 20);
  o += feed_forward(s, nn_human_KIN_RSK_group_2, 13, 10);
  o += feed_forward(s, nn_human_KIN_RSK_group_3, 13, 15);
  o += feed_forward(s, nn_human_KIN_RSK_group_4, 13, 10);
  o += feed_forward(s, nn_human_KIN_RSK_group_5, 5, 15);
  o += feed_forward(s, nn_human_KIN_RSK_group_6, 11, 0);
  o += feed_forward(s, nn_human_KIN_RSK_group_7, 7, 2);
  o += feed_forward(s, nn_human_KIN_RSK_group_8, 7, 4);
  o += feed_forward(s, nn_human_KIN_RSK_group_9, 7, 6);
  o += feed_forward(s, nn_human_KIN_RSK_group_10, 7, 20);
  o += feed_forward(s, nn_human_KIN_RSK_group_11, 13, 10);
  o /= 11;
  o = 0.024199722457706+(0.10878128746544-0.024199722457706)/(1+exp(99.9999*(0.505252-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tRSK_group\t%.6f\t%.6f\n", o, 0.04);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_SGK_group_1, 13, 0);
  o += feed_forward(s, nn_human_KIN_SGK_group_2, 11, 0);
  o += feed_forward(s, nn_human_KIN_SGK_group_3, 7, 6);
  o += feed_forward(s, nn_human_KIN_SGK_group_4, 11, 15);
  o += feed_forward(s, nn_human_KIN_SGK_group_5, 7, 10);
  o += feed_forward(s, nn_human_KIN_SGK_group_6, 9, 2);
  o += feed_forward(s, nn_human_KIN_SGK_group_7, 9, 0);
  o += feed_forward(s, nn_human_KIN_SGK_group_8, 7, 20);
  o += feed_forward(s, nn_human_KIN_SGK_group_9, 9, 10);
  o += feed_forward(s, nn_human_KIN_SGK_group_10, 11, 2);
  o += feed_forward(s, nn_human_KIN_SGK_group_11, 9, 0);
  o += feed_forward(s, nn_human_KIN_SGK_group_12, 9, 4);
  o /= 12;
  o = 0.0131074375261506+(0.123704261566443-0.0131074375261506)/(1+exp(100*(0.505722-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tSGK_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_Syk_group_1, 9, 20);
  o += feed_forward(s, nn_human_KIN_Syk_group_2, 9, 2);
  o += feed_forward(s, nn_human_KIN_Syk_group_3, 5, 4);
  o += feed_forward(s, nn_human_KIN_Syk_group_4, 13, 4);
  o += feed_forward(s, nn_human_KIN_Syk_group_5, 5, 15);
  o += feed_forward(s, nn_human_KIN_Syk_group_6, 13, 6);
  o += feed_forward(s, nn_human_KIN_Syk_group_7, 11, 4);
  o += feed_forward(s, nn_human_KIN_Syk_group_8, 7, 6);
  o += feed_forward(s, nn_human_KIN_Syk_group_9, 13, 20);
  o += feed_forward(s, nn_human_KIN_Syk_group_10, 11, 20);
  o += feed_forward(s, nn_human_KIN_Syk_group_11, 9, 10);
  o += feed_forward(s, nn_human_KIN_Syk_group_12, 7, 20);
  o /= 12;
  o = 0.0159425432882374+(0.0792141690740016-0.0159425432882374)/(1+exp(92.3993*(0.497051-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tSyk_group\t%.6f\t%.6f\n", o, 0.0282842712474619);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_TTK_1, 9, 10);
  o += feed_forward(s, nn_human_KIN_TTK_2, 5, 10);
  o += feed_forward(s, nn_human_KIN_TTK_3, 5, 2);
  o += feed_forward(s, nn_human_KIN_TTK_4, 9, 20);
  o += feed_forward(s, nn_human_KIN_TTK_5, 5, 6);
  o += feed_forward(s, nn_human_KIN_TTK_6, 5, 10);
  o += feed_forward(s, nn_human_KIN_TTK_7, 9, 20);
  o += feed_forward(s, nn_human_KIN_TTK_8, 7, 6);
  o += feed_forward(s, nn_human_KIN_TTK_9, 5, 10);
  o += feed_forward(s, nn_human_KIN_TTK_10, 9, 20);
  o += feed_forward(s, nn_human_KIN_TTK_11, 9, 10);
  o += feed_forward(s, nn_human_KIN_TTK_12, 5, 10);
  o /= 12;
  o = 7.59970842354081e-09+(0.0601937269372694-7.59970842354081e-09)/(1+exp(100*(0.493895-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tTTK\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_Tec_group_1, 5, 0);
  o += feed_forward(s, nn_human_KIN_Tec_group_2, 13, 2);
  o += feed_forward(s, nn_human_KIN_Tec_group_3, 7, 10);
  o += feed_forward(s, nn_human_KIN_Tec_group_4, 7, 2);
  o += feed_forward(s, nn_human_KIN_Tec_group_5, 13, 2);
  o += feed_forward(s, nn_human_KIN_Tec_group_6, 13, 0);
  o += feed_forward(s, nn_human_KIN_Tec_group_7, 5, 2);
  o += feed_forward(s, nn_human_KIN_Tec_group_8, 13, 2);
  o += feed_forward(s, nn_human_KIN_Tec_group_9, 7, 10);
  o += feed_forward(s, nn_human_KIN_Tec_group_10, 7, 10);
  o += feed_forward(s, nn_human_KIN_Tec_group_11, 7, 15);
  o /= 11;
  o = 7.9720878967667e-05+(0.0801631521630341-7.9720878967667e-05)/(1+exp(31.3671*(0.484583-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tTec_group\t%.6f\t%.6f\n", o, 0.0447213595499958);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_Trk_group_1, 9, 10);
  o += feed_forward(s, nn_human_KIN_Trk_group_2, 13, 15);
  o += feed_forward(s, nn_human_KIN_Trk_group_3, 9, 0);
  o += feed_forward(s, nn_human_KIN_Trk_group_4, 13, 4);
  o += feed_forward(s, nn_human_KIN_Trk_group_5, 5, 6);
  o += feed_forward(s, nn_human_KIN_Trk_group_6, 5, 2);
  o += feed_forward(s, nn_human_KIN_Trk_group_7, 7, 15);
  o += feed_forward(s, nn_human_KIN_Trk_group_8, 13, 4);
  o += feed_forward(s, nn_human_KIN_Trk_group_9, 5, 6);
  o += feed_forward(s, nn_human_KIN_Trk_group_10, 13, 15);
  o += feed_forward(s, nn_human_KIN_Trk_group_11, 13, 10);
  o += feed_forward(s, nn_human_KIN_Trk_group_12, 7, 15);
  o /= 12;
  o = 2.19586271302407e-08+(0.103421429055034-2.19586271302407e-08)/(1+exp(99.9999*(0.494324-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tTrk_group\t%.6f\t%.6f\n", o, 0.0346410161513775);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_KIN_Tyk2_1, 5, 10);
  o += feed_forward(s, nn_human_KIN_Tyk2_2, 5, 20);
  o += feed_forward(s, nn_human_KIN_Tyk2_3, 7, 10);
  o += feed_forward(s, nn_human_KIN_Tyk2_4, 7, 2);
  o += feed_forward(s, nn_human_KIN_Tyk2_5, 7, 2);
  o += feed_forward(s, nn_human_KIN_Tyk2_6, 5, 0);
  o += feed_forward(s, nn_human_KIN_Tyk2_7, 7, 15);
  o += feed_forward(s, nn_human_KIN_Tyk2_8, 7, 2);
  o /= 8;
  o = 0.00746772424267994+(0.0475933520397151-0.00746772424267994)/(1+exp(100*(0.502659-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tKIN\tTyk2\t%.6f\t%.6f\n", o, 0.02);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_NT1_group_1, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT1_group_2, 5, 0);
  o += feed_forward(s, nn_human_PTP_NT1_group_3, 7, 2);
  o += feed_forward(s, nn_human_PTP_NT1_group_4, 5, 2);
  o += feed_forward(s, nn_human_PTP_NT1_group_5, 9, 2);
  o += feed_forward(s, nn_human_PTP_NT1_group_6, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT1_group_7, 5, 2);
  o += feed_forward(s, nn_human_PTP_NT1_group_8, 9, 0);
  o += feed_forward(s, nn_human_PTP_NT1_group_9, 7, 2);
  o += feed_forward(s, nn_human_PTP_NT1_group_10, 9, 0);
  o /= 10;
  o = 0.0109319410282579+(0.20655123892423-0.0109319410282579)/(1+exp(10.6207*(0.573087-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tNT1_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_NT2_group_1, 7, 2);
  o += feed_forward(s, nn_human_PTP_NT2_group_2, 5, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_3, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_4, 7, 2);
  o += feed_forward(s, nn_human_PTP_NT2_group_5, 9, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_6, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_7, 9, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_8, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_9, 9, 0);
  o += feed_forward(s, nn_human_PTP_NT2_group_10, 5, 0);
  o /= 10;
  o = 0.0261848455144826+(0.0791496430680587-0.0261848455144826)/(1+exp(30.9822*(0.493289-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tNT2_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_NT4_group_1, 5, 0);
  o += feed_forward(s, nn_human_PTP_NT4_group_2, 7, 2);
  o += feed_forward(s, nn_human_PTP_NT4_group_3, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT4_group_4, 9, 2);
  o += feed_forward(s, nn_human_PTP_NT4_group_5, 9, 0);
  o += feed_forward(s, nn_human_PTP_NT4_group_6, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT4_group_7, 9, 0);
  o += feed_forward(s, nn_human_PTP_NT4_group_8, 7, 0);
  o += feed_forward(s, nn_human_PTP_NT4_group_9, 9, 2);
  o += feed_forward(s, nn_human_PTP_NT4_group_10, 7, 2);
  o /= 10;
  o = 0.00686121310803699+(0.181698699565983-0.00686121310803699)/(1+exp(9.69547*(0.520943-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tNT4_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_PTPN13_1, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN13_2, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN13_3, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN13_4, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN13_5, 7, 2);
  o /= 5;
  o = 0.0148296092069492+(0.0700001627153945-0.0148296092069492)/(1+exp(50.0458*(0.51502-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tPTPN13\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_PTPN23_1, 9, 0);
  o += feed_forward(s, nn_human_PTP_PTPN23_2, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN23_3, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN23_4, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN23_5, 7, 0);
  o += feed_forward(s, nn_human_PTP_PTPN23_6, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN23_7, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN23_8, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN23_9, 7, 0);
  o /= 9;
  o = 9.86843492811949e-07+(0.103521878335112-9.86843492811949e-07)/(1+exp(11.9027*(0.465505-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tPTPN23\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_PTPN3_1, 5, 2);
  o += feed_forward(s, nn_human_PTP_PTPN3_2, 7, 0);
  o += feed_forward(s, nn_human_PTP_PTPN3_3, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN3_4, 5, 2);
  o += feed_forward(s, nn_human_PTP_PTPN3_5, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN3_6, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN3_7, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN3_8, 9, 0);
  o += feed_forward(s, nn_human_PTP_PTPN3_9, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN3_10, 5, 2);
  o /= 10;
  o = 0.0125998408600296+(0.0612606654447649-0.0125998408600296)/(1+exp(47.2508*(0.500298-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tPTPN3\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_PTPN4_1, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN4_2, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN4_3, 5, 2);
  o += feed_forward(s, nn_human_PTP_PTPN4_4, 7, 0);
  o += feed_forward(s, nn_human_PTP_PTPN4_5, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN4_6, 5, 2);
  o += feed_forward(s, nn_human_PTP_PTPN4_7, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN4_8, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN4_9, 9, 2);
  o /= 9;
  o = 0.00649515073940185+(0.134716111460297-0.00649515073940185)/(1+exp(19.3334*(0.510239-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tPTPN4\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_PTPN9_1, 7, 2);
  o += feed_forward(s, nn_human_PTP_PTPN9_2, 5, 2);
  o += feed_forward(s, nn_human_PTP_PTPN9_3, 9, 2);
  o += feed_forward(s, nn_human_PTP_PTPN9_4, 7, 0);
  o += feed_forward(s, nn_human_PTP_PTPN9_5, 5, 2);
  o += feed_forward(s, nn_human_PTP_PTPN9_6, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN9_7, 5, 0);
  o += feed_forward(s, nn_human_PTP_PTPN9_8, 9, 2);
  o /= 8;
  o = 0.00706204772090403+(0.1426784802999-0.00706204772090403)/(1+exp(10.5751*(0.561225-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tPTPN9\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_R2A_group_1, 7, 2);
  o += feed_forward(s, nn_human_PTP_R2A_group_2, 5, 2);
  o += feed_forward(s, nn_human_PTP_R2A_group_3, 7, 0);
  o += feed_forward(s, nn_human_PTP_R2A_group_4, 9, 0);
  o += feed_forward(s, nn_human_PTP_R2A_group_5, 9, 2);
  o += feed_forward(s, nn_human_PTP_R2A_group_6, 9, 0);
  o += feed_forward(s, nn_human_PTP_R2A_group_7, 5, 2);
  o += feed_forward(s, nn_human_PTP_R2A_group_8, 9, 2);
  o += feed_forward(s, nn_human_PTP_R2A_group_9, 5, 2);
  o += feed_forward(s, nn_human_PTP_R2A_group_10, 9, 0);
  o /= 10;
  o = 0.0270141944363556+(0.167195764855025-0.0270141944363556)/(1+exp(31.0353*(0.514323-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tR2A_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_R3_group_1, 5, 0);
  o += feed_forward(s, nn_human_PTP_R3_group_2, 7, 2);
  o += feed_forward(s, nn_human_PTP_R3_group_3, 5, 2);
  o += feed_forward(s, nn_human_PTP_R3_group_4, 5, 0);
  o += feed_forward(s, nn_human_PTP_R3_group_5, 7, 2);
  o += feed_forward(s, nn_human_PTP_R3_group_6, 7, 0);
  o += feed_forward(s, nn_human_PTP_R3_group_7, 5, 2);
  o += feed_forward(s, nn_human_PTP_R3_group_8, 5, 0);
  o += feed_forward(s, nn_human_PTP_R3_group_9, 7, 2);
  o += feed_forward(s, nn_human_PTP_R3_group_10, 5, 0);
  o += feed_forward(s, nn_human_PTP_R3_group_11, 9, 2);
  o += feed_forward(s, nn_human_PTP_R3_group_12, 7, 0);
  o /= 12;
  o = 0.0112805237146122+(0.422575629185381-0.0112805237146122)/(1+exp(9.89258*(0.526406-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tR3_group\t%.6f\t%.6f\n", o, 0.063887656499994);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTP_R4_group_1, 7, 0);
  o += feed_forward(s, nn_human_PTP_R4_group_2, 9, 2);
  o += feed_forward(s, nn_human_PTP_R4_group_3, 9, 2);
  o += feed_forward(s, nn_human_PTP_R4_group_4, 5, 0);
  o += feed_forward(s, nn_human_PTP_R4_group_5, 7, 2);
  o += feed_forward(s, nn_human_PTP_R4_group_6, 5, 2);
  o += feed_forward(s, nn_human_PTP_R4_group_7, 7, 2);
  o += feed_forward(s, nn_human_PTP_R4_group_8, 7, 0);
  o += feed_forward(s, nn_human_PTP_R4_group_9, 7, 2);
  o += feed_forward(s, nn_human_PTP_R4_group_10, 9, 0);
  o += feed_forward(s, nn_human_PTP_R4_group_11, 5, 2);
  o /= 11;
  o = 0.0135180705374386+(0.149184049063479-0.0135180705374386)/(1+exp(12.2294*(0.516195-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTP\tR4_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Abl_group_1, 11, 15);
  o += feed_forward(s, nn_human_SH2_Abl_group_2, 9, 2);
  o += feed_forward(s, nn_human_SH2_Abl_group_3, 13, 20);
  o += feed_forward(s, nn_human_SH2_Abl_group_4, 13, 4);
  o += feed_forward(s, nn_human_SH2_Abl_group_5, 7, 0);
  o += feed_forward(s, nn_human_SH2_Abl_group_6, 13, 4);
  o += feed_forward(s, nn_human_SH2_Abl_group_7, 13, 6);
  o += feed_forward(s, nn_human_SH2_Abl_group_8, 9, 20);
  o += feed_forward(s, nn_human_SH2_Abl_group_9, 7, 4);
  o += feed_forward(s, nn_human_SH2_Abl_group_10, 13, 10);
  o += feed_forward(s, nn_human_SH2_Abl_group_11, 13, 4);
  o += feed_forward(s, nn_human_SH2_Abl_group_12, 13, 4);
  o /= 12;
  o = 0.0126971577412434+(0.162721014323648-0.0126971577412434)/(1+exp(12.489*(0.521011-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tAbl_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_BCAR3_1, 11, 4);
  o += feed_forward(s, nn_human_SH2_BCAR3_2, 13, 4);
  o += feed_forward(s, nn_human_SH2_BCAR3_3, 13, 0);
  o += feed_forward(s, nn_human_SH2_BCAR3_4, 11, 2);
  o += feed_forward(s, nn_human_SH2_BCAR3_5, 13, 4);
  o += feed_forward(s, nn_human_SH2_BCAR3_6, 11, 20);
  o += feed_forward(s, nn_human_SH2_BCAR3_7, 7, 15);
  o += feed_forward(s, nn_human_SH2_BCAR3_8, 7, 15);
  o += feed_forward(s, nn_human_SH2_BCAR3_9, 7, 15);
  o += feed_forward(s, nn_human_SH2_BCAR3_10, 11, 20);
  o += feed_forward(s, nn_human_SH2_BCAR3_11, 9, 10);
  o += feed_forward(s, nn_human_SH2_BCAR3_12, 7, 20);
  o /= 12;
  o = 0.0141104597534641+(0.0733241515406911-0.0141104597534641)/(1+exp(26.5905*(0.486638-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tBCAR3\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_BRDG1_1, 9, 0);
  o += feed_forward(s, nn_human_SH2_BRDG1_2, 9, 20);
  o += feed_forward(s, nn_human_SH2_BRDG1_3, 13, 10);
  o += feed_forward(s, nn_human_SH2_BRDG1_4, 11, 20);
  o += feed_forward(s, nn_human_SH2_BRDG1_5, 11, 20);
  o += feed_forward(s, nn_human_SH2_BRDG1_6, 9, 20);
  o += feed_forward(s, nn_human_SH2_BRDG1_7, 11, 6);
  o += feed_forward(s, nn_human_SH2_BRDG1_8, 9, 6);
  o += feed_forward(s, nn_human_SH2_BRDG1_9, 13, 10);
  o += feed_forward(s, nn_human_SH2_BRDG1_10, 11, 0);
  o += feed_forward(s, nn_human_SH2_BRDG1_11, 11, 6);
  o += feed_forward(s, nn_human_SH2_BRDG1_12, 9, 2);
  o /= 12;
  o = 0.00394723879722755+(0.211191523995335-0.00394723879722755)/(1+exp(23.7797*(0.490839-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tBRDG1\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Brk_1, 7, 0);
  o += feed_forward(s, nn_human_SH2_Brk_2, 5, 20);
  o += feed_forward(s, nn_human_SH2_Brk_3, 7, 2);
  o += feed_forward(s, nn_human_SH2_Brk_4, 13, 15);
  o += feed_forward(s, nn_human_SH2_Brk_5, 5, 20);
  o += feed_forward(s, nn_human_SH2_Brk_6, 7, 4);
  o += feed_forward(s, nn_human_SH2_Brk_7, 9, 6);
  o += feed_forward(s, nn_human_SH2_Brk_8, 7, 0);
  o += feed_forward(s, nn_human_SH2_Brk_9, 5, 10);
  o += feed_forward(s, nn_human_SH2_Brk_10, 7, 20);
  o += feed_forward(s, nn_human_SH2_Brk_11, 13, 0);
  o += feed_forward(s, nn_human_SH2_Brk_12, 5, 20);
  o /= 12;
  o = 0.00350879038993763+(0.202952029520295-0.00350879038993763)/(1+exp(17.5139*(0.535357-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tBrk\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_CBL_group_1, 13, 15);
  o += feed_forward(s, nn_human_SH2_CBL_group_2, 11, 15);
  o += feed_forward(s, nn_human_SH2_CBL_group_3, 11, 15);
  o += feed_forward(s, nn_human_SH2_CBL_group_4, 5, 4);
  o += feed_forward(s, nn_human_SH2_CBL_group_5, 7, 6);
  o += feed_forward(s, nn_human_SH2_CBL_group_6, 9, 20);
  o += feed_forward(s, nn_human_SH2_CBL_group_7, 11, 20);
  o += feed_forward(s, nn_human_SH2_CBL_group_8, 11, 15);
  o += feed_forward(s, nn_human_SH2_CBL_group_9, 13, 15);
  o += feed_forward(s, nn_human_SH2_CBL_group_10, 9, 0);
  o += feed_forward(s, nn_human_SH2_CBL_group_11, 13, 4);
  o += feed_forward(s, nn_human_SH2_CBL_group_12, 5, 6);
  o /= 12;
  o = 0.0175830289649631+(0.148900717575677-0.0175830289649631)/(1+exp(21.1237*(0.488658-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tCBL_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_CRK_1, 13, 4);
  o += feed_forward(s, nn_human_SH2_CRK_2, 7, 2);
  o += feed_forward(s, nn_human_SH2_CRK_3, 7, 6);
  o += feed_forward(s, nn_human_SH2_CRK_4, 13, 20);
  o += feed_forward(s, nn_human_SH2_CRK_5, 11, 0);
  o += feed_forward(s, nn_human_SH2_CRK_6, 7, 6);
  o += feed_forward(s, nn_human_SH2_CRK_7, 13, 6);
  o += feed_forward(s, nn_human_SH2_CRK_8, 13, 10);
  o += feed_forward(s, nn_human_SH2_CRK_9, 7, 6);
  o += feed_forward(s, nn_human_SH2_CRK_10, 13, 6);
  o += feed_forward(s, nn_human_SH2_CRK_11, 9, 10);
  o /= 11;
  o = 0.00385099479719505+(0.239646887770349-0.00385099479719505)/(1+exp(15.7559*(0.521387-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tCRK\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_CRKL_1, 11, 2);
  o += feed_forward(s, nn_human_SH2_CRKL_2, 9, 0);
  o += feed_forward(s, nn_human_SH2_CRKL_3, 7, 10);
  o += feed_forward(s, nn_human_SH2_CRKL_4, 13, 15);
  o += feed_forward(s, nn_human_SH2_CRKL_5, 11, 15);
  o += feed_forward(s, nn_human_SH2_CRKL_6, 9, 6);
  o += feed_forward(s, nn_human_SH2_CRKL_7, 9, 2);
  o += feed_forward(s, nn_human_SH2_CRKL_8, 7, 6);
  o += feed_forward(s, nn_human_SH2_CRKL_9, 11, 15);
  o += feed_forward(s, nn_human_SH2_CRKL_10, 9, 4);
  o += feed_forward(s, nn_human_SH2_CRKL_11, 11, 2);
  o /= 11;
  o = 0.00343164763876544+(0.302105612477238-0.00343164763876544)/(1+exp(17.6546*(0.515251-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tCRKL\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_CTK_1, 11, 2);
  o += feed_forward(s, nn_human_SH2_CTK_2, 13, 10);
  o += feed_forward(s, nn_human_SH2_CTK_3, 13, 2);
  o += feed_forward(s, nn_human_SH2_CTK_4, 5, 10);
  o += feed_forward(s, nn_human_SH2_CTK_5, 9, 4);
  o += feed_forward(s, nn_human_SH2_CTK_6, 7, 2);
  o += feed_forward(s, nn_human_SH2_CTK_7, 11, 2);
  o += feed_forward(s, nn_human_SH2_CTK_8, 5, 6);
  o += feed_forward(s, nn_human_SH2_CTK_9, 7, 6);
  o += feed_forward(s, nn_human_SH2_CTK_10, 5, 2);
  o += feed_forward(s, nn_human_SH2_CTK_11, 5, 15);
  o += feed_forward(s, nn_human_SH2_CTK_12, 13, 15);
  o /= 12;
  o = 0.00960347439939094+(0.092460889668254-0.00960347439939094)/(1+exp(99.891*(0.472228-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tCTK\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Csk_1, 11, 6);
  o += feed_forward(s, nn_human_SH2_Csk_2, 13, 10);
  o += feed_forward(s, nn_human_SH2_Csk_3, 11, 10);
  o += feed_forward(s, nn_human_SH2_Csk_4, 11, 15);
  o += feed_forward(s, nn_human_SH2_Csk_5, 13, 0);
  o += feed_forward(s, nn_human_SH2_Csk_6, 9, 20);
  o += feed_forward(s, nn_human_SH2_Csk_7, 13, 15);
  o += feed_forward(s, nn_human_SH2_Csk_8, 11, 15);
  o += feed_forward(s, nn_human_SH2_Csk_9, 9, 20);
  o += feed_forward(s, nn_human_SH2_Csk_10, 11, 0);
  o += feed_forward(s, nn_human_SH2_Csk_11, 13, 20);
  o += feed_forward(s, nn_human_SH2_Csk_12, 11, 20);
  o /= 12;
  o = 0.00787233041956804+(0.175839822461384-0.00787233041956804)/(1+exp(28.4319*(0.503179-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tCsk\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_DAPP1_1, 5, 15);
  o += feed_forward(s, nn_human_SH2_DAPP1_2, 5, 20);
  o += feed_forward(s, nn_human_SH2_DAPP1_3, 11, 6);
  o += feed_forward(s, nn_human_SH2_DAPP1_4, 5, 0);
  o += feed_forward(s, nn_human_SH2_DAPP1_5, 5, 20);
  o += feed_forward(s, nn_human_SH2_DAPP1_6, 5, 6);
  o += feed_forward(s, nn_human_SH2_DAPP1_7, 11, 20);
  o += feed_forward(s, nn_human_SH2_DAPP1_8, 7, 15);
  o += feed_forward(s, nn_human_SH2_DAPP1_9, 7, 2);
  o += feed_forward(s, nn_human_SH2_DAPP1_10, 13, 4);
  o += feed_forward(s, nn_human_SH2_DAPP1_11, 5, 15);
  o += feed_forward(s, nn_human_SH2_DAPP1_12, 9, 10);
  o /= 12;
  o = 0.0171865116316951+(0.0826855123674912-0.0171865116316951)/(1+exp(100*(0.501927-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tDAPP1\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_FRK_1, 5, 0);
  o += feed_forward(s, nn_human_SH2_FRK_2, 9, 2);
  o += feed_forward(s, nn_human_SH2_FRK_3, 11, 4);
  o += feed_forward(s, nn_human_SH2_FRK_4, 9, 20);
  o += feed_forward(s, nn_human_SH2_FRK_5, 7, 2);
  o += feed_forward(s, nn_human_SH2_FRK_6, 11, 15);
  o += feed_forward(s, nn_human_SH2_FRK_7, 11, 4);
  o += feed_forward(s, nn_human_SH2_FRK_8, 11, 2);
  o += feed_forward(s, nn_human_SH2_FRK_9, 11, 15);
  o += feed_forward(s, nn_human_SH2_FRK_10, 5, 6);
  o += feed_forward(s, nn_human_SH2_FRK_11, 11, 15);
  o /= 11;
  o = 0.010672433572007+(0.120938646099411-0.010672433572007)/(1+exp(88.4599*(0.496526-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tFRK\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Fes_group_1, 9, 6);
  o += feed_forward(s, nn_human_SH2_Fes_group_2, 13, 6);
  o += feed_forward(s, nn_human_SH2_Fes_group_3, 9, 6);
  o += feed_forward(s, nn_human_SH2_Fes_group_4, 7, 10);
  o += feed_forward(s, nn_human_SH2_Fes_group_5, 9, 20);
  o += feed_forward(s, nn_human_SH2_Fes_group_6, 7, 20);
  o += feed_forward(s, nn_human_SH2_Fes_group_7, 9, 20);
  o += feed_forward(s, nn_human_SH2_Fes_group_8, 9, 15);
  o += feed_forward(s, nn_human_SH2_Fes_group_9, 13, 2);
  o += feed_forward(s, nn_human_SH2_Fes_group_10, 7, 20);
  o += feed_forward(s, nn_human_SH2_Fes_group_11, 9, 15);
  o += feed_forward(s, nn_human_SH2_Fes_group_12, 9, 15);
  o /= 12;
  o = 0.00696923516146076+(0.294228088765004-0.00696923516146076)/(1+exp(38.2392*(0.486657-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tFes_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_GRB2_group_1, 7, 2);
  o += feed_forward(s, nn_human_SH2_GRB2_group_2, 5, 10);
  o += feed_forward(s, nn_human_SH2_GRB2_group_3, 11, 6);
  o += feed_forward(s, nn_human_SH2_GRB2_group_4, 11, 10);
  o += feed_forward(s, nn_human_SH2_GRB2_group_5, 11, 6);
  o += feed_forward(s, nn_human_SH2_GRB2_group_6, 9, 0);
  o += feed_forward(s, nn_human_SH2_GRB2_group_7, 11, 2);
  o += feed_forward(s, nn_human_SH2_GRB2_group_8, 11, 6);
  o += feed_forward(s, nn_human_SH2_GRB2_group_9, 13, 20);
  o += feed_forward(s, nn_human_SH2_GRB2_group_10, 5, 10);
  o += feed_forward(s, nn_human_SH2_GRB2_group_11, 11, 2);
  o += feed_forward(s, nn_human_SH2_GRB2_group_12, 9, 20);
  o /= 12;
  o = 1.31995534340683e-05+(0.5773684663168-1.31995534340683e-05)/(1+exp(22.6919*(0.511021-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tGRB2_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_GRB_group_1, 13, 20);
  o += feed_forward(s, nn_human_SH2_GRB_group_2, 13, 10);
  o += feed_forward(s, nn_human_SH2_GRB_group_3, 11, 15);
  o += feed_forward(s, nn_human_SH2_GRB_group_4, 13, 20);
  o += feed_forward(s, nn_human_SH2_GRB_group_5, 7, 6);
  o += feed_forward(s, nn_human_SH2_GRB_group_6, 5, 4);
  o += feed_forward(s, nn_human_SH2_GRB_group_7, 7, 6);
  o += feed_forward(s, nn_human_SH2_GRB_group_8, 5, 2);
  o += feed_forward(s, nn_human_SH2_GRB_group_9, 5, 4);
  o += feed_forward(s, nn_human_SH2_GRB_group_10, 7, 20);
  o += feed_forward(s, nn_human_SH2_GRB_group_11, 5, 20);
  o += feed_forward(s, nn_human_SH2_GRB_group_12, 7, 15);
  o /= 12;
  o = 0.026396349969728+(0.14927182749983-0.026396349969728)/(1+exp(44.0508*(0.50202-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tGRB_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_1, 13, 15);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_2, 13, 6);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_3, 9, 10);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_4, 9, 2);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_5, 9, 10);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_6, 11, 4);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_7, 7, 10);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_8, 9, 4);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_9, 9, 10);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_10, 9, 10);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_11, 9, 15);
  o += feed_forward(s, nn_human_SH2_HSH2D_SH2D2A_group_12, 5, 15);
  o /= 12;
  o = 0.0122600207314021+(0.143151040791059-0.0122600207314021)/(1+exp(99.9887*(0.491558-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tHSH2D_SH2D2A_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_INPPL1_1, 9, 15);
  o += feed_forward(s, nn_human_SH2_INPPL1_2, 7, 2);
  o += feed_forward(s, nn_human_SH2_INPPL1_3, 9, 10);
  o += feed_forward(s, nn_human_SH2_INPPL1_4, 11, 4);
  o += feed_forward(s, nn_human_SH2_INPPL1_5, 9, 15);
  o += feed_forward(s, nn_human_SH2_INPPL1_6, 13, 2);
  o += feed_forward(s, nn_human_SH2_INPPL1_7, 7, 10);
  o += feed_forward(s, nn_human_SH2_INPPL1_8, 11, 6);
  o += feed_forward(s, nn_human_SH2_INPPL1_9, 9, 15);
  o += feed_forward(s, nn_human_SH2_INPPL1_10, 7, 2);
  o += feed_forward(s, nn_human_SH2_INPPL1_11, 5, 15);
  o += feed_forward(s, nn_human_SH2_INPPL1_12, 9, 0);
  o /= 12;
  o = 9.65859830830319e-06+(0.0863948966789007-9.65859830830319e-06)/(1+exp(8.44492*(0.449336-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tINPPL1\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_1, 13, 15);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_2, 7, 2);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_3, 9, 20);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_4, 11, 20);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_5, 5, 15);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_6, 9, 20);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_7, 11, 15);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_8, 9, 20);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_9, 9, 10);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_10, 7, 2);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_11, 5, 20);
  o += feed_forward(s, nn_human_SH2_MIST_LCP2_BLNK_SH2D6_group_12, 7, 10);
  o /= 12;
  o = 0.0282117919528736+(0.251414507564104-0.0282117919528736)/(1+exp(94.661*(0.496338-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tMIST_LCP2_BLNK_SH2D6_group\t%.6f\t%.6f\n", o, 0.0571428571428571);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_NCK_group_1, 7, 20);
  o += feed_forward(s, nn_human_SH2_NCK_group_2, 7, 6);
  o += feed_forward(s, nn_human_SH2_NCK_group_3, 11, 15);
  o += feed_forward(s, nn_human_SH2_NCK_group_4, 9, 4);
  o += feed_forward(s, nn_human_SH2_NCK_group_5, 7, 0);
  o += feed_forward(s, nn_human_SH2_NCK_group_6, 11, 20);
  o += feed_forward(s, nn_human_SH2_NCK_group_7, 7, 2);
  o += feed_forward(s, nn_human_SH2_NCK_group_8, 7, 6);
  o += feed_forward(s, nn_human_SH2_NCK_group_9, 11, 20);
  o += feed_forward(s, nn_human_SH2_NCK_group_10, 9, 10);
  o += feed_forward(s, nn_human_SH2_NCK_group_11, 7, 15);
  o += feed_forward(s, nn_human_SH2_NCK_group_12, 11, 6);
  o /= 12;
  o = 1.01751942109341e-06+(0.163972208551221-1.01751942109341e-06)/(1+exp(16.2913*(0.469711-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tNCK_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_1, 7, 20);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_2, 7, 15);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_3, 9, 2);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_4, 13, 10);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_5, 7, 10);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_6, 7, 4);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_7, 11, 2);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_8, 9, 15);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_9, 7, 20);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_10, 9, 6);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_11, 11, 10);
  o += feed_forward(s, nn_human_SH2_PIK3R1_1_12, 7, 0);
  o /= 12;
  o = 0.00625004164957589+(0.396779402454771-0.00625004164957589)/(1+exp(19.0055*(0.544806-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tPIK3R1_1\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_1, 9, 0);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_2, 7, 10);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_3, 9, 6);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_4, 9, 20);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_5, 7, 15);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_6, 13, 6);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_7, 11, 4);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_8, 13, 20);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_9, 7, 10);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_10, 7, 2);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_11, 11, 15);
  o += feed_forward(s, nn_human_SH2_PIK3R3_1_PIK3R2_1_group_12, 7, 0);
  o /= 12;
  o = 0.00640690060102454+(0.544627563507808-0.00640690060102454)/(1+exp(26.0913*(0.524418-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tPIK3R3_1_PIK3R2_1_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_PTPN_group_1, 7, 4);
  o += feed_forward(s, nn_human_SH2_PTPN_group_2, 7, 15);
  o += feed_forward(s, nn_human_SH2_PTPN_group_3, 13, 0);
  o += feed_forward(s, nn_human_SH2_PTPN_group_4, 11, 15);
  o += feed_forward(s, nn_human_SH2_PTPN_group_5, 5, 0);
  o += feed_forward(s, nn_human_SH2_PTPN_group_6, 11, 20);
  o += feed_forward(s, nn_human_SH2_PTPN_group_7, 7, 20);
  o += feed_forward(s, nn_human_SH2_PTPN_group_8, 7, 4);
  o += feed_forward(s, nn_human_SH2_PTPN_group_9, 9, 20);
  o += feed_forward(s, nn_human_SH2_PTPN_group_10, 11, 15);
  o += feed_forward(s, nn_human_SH2_PTPN_group_11, 13, 0);
  o += feed_forward(s, nn_human_SH2_PTPN_group_12, 9, 4);
  o /= 12;
  o = 0.00886336766617263+(0.283547528076879-0.00886336766617263)/(1+exp(14.2606*(0.521992-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tPTPN_group\t%.6f\t%.6f\n", o, 0.0571428571428571);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_RASA_group_1, 9, 10);
  o += feed_forward(s, nn_human_SH2_RASA_group_2, 11, 2);
  o += feed_forward(s, nn_human_SH2_RASA_group_3, 11, 4);
  o += feed_forward(s, nn_human_SH2_RASA_group_4, 11, 2);
  o += feed_forward(s, nn_human_SH2_RASA_group_5, 7, 2);
  o += feed_forward(s, nn_human_SH2_RASA_group_6, 11, 15);
  o += feed_forward(s, nn_human_SH2_RASA_group_7, 7, 15);
  o += feed_forward(s, nn_human_SH2_RASA_group_8, 11, 2);
  o += feed_forward(s, nn_human_SH2_RASA_group_9, 7, 15);
  o += feed_forward(s, nn_human_SH2_RASA_group_10, 13, 4);
  o += feed_forward(s, nn_human_SH2_RASA_group_11, 7, 2);
  o += feed_forward(s, nn_human_SH2_RASA_group_12, 7, 20);
  o /= 12;
  o = 0.0123637140633305+(0.201572324786627-0.0123637140633305)/(1+exp(61.1431*(0.503898-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tRASA_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SH2B_group_1, 7, 6);
  o += feed_forward(s, nn_human_SH2_SH2B_group_2, 11, 15);
  o += feed_forward(s, nn_human_SH2_SH2B_group_3, 9, 4);
  o += feed_forward(s, nn_human_SH2_SH2B_group_4, 7, 10);
  o += feed_forward(s, nn_human_SH2_SH2B_group_5, 11, 15);
  o += feed_forward(s, nn_human_SH2_SH2B_group_6, 11, 10);
  o += feed_forward(s, nn_human_SH2_SH2B_group_7, 7, 0);
  o += feed_forward(s, nn_human_SH2_SH2B_group_8, 7, 15);
  o += feed_forward(s, nn_human_SH2_SH2B_group_9, 7, 6);
  o += feed_forward(s, nn_human_SH2_SH2B_group_10, 9, 6);
  o += feed_forward(s, nn_human_SH2_SH2B_group_11, 13, 2);
  o += feed_forward(s, nn_human_SH2_SH2B_group_12, 13, 20);
  o /= 12;
  o = 0.0202426008107616+(0.194384506485122-0.0202426008107616)/(1+exp(21.007*(0.526091-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSH2B_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_1, 11, 6);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_2, 7, 15);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_3, 11, 15);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_4, 13, 4);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_5, 7, 0);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_6, 7, 15);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_7, 9, 6);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_8, 11, 15);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_9, 13, 10);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_10, 11, 0);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_11, 13, 2);
  o += feed_forward(s, nn_human_SH2_SH2D1A_SH2D1B_group_12, 9, 15);
  o /= 12;
  o = 0.000391769586359461+(0.19445959903907-0.000391769586359461)/(1+exp(17.4598*(0.50565-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSH2D1A_SH2D1B_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SH2D3C_1, 7, 20);
  o += feed_forward(s, nn_human_SH2_SH2D3C_2, 11, 20);
  o += feed_forward(s, nn_human_SH2_SH2D3C_3, 11, 2);
  o += feed_forward(s, nn_human_SH2_SH2D3C_4, 9, 10);
  o += feed_forward(s, nn_human_SH2_SH2D3C_5, 7, 15);
  o += feed_forward(s, nn_human_SH2_SH2D3C_6, 7, 20);
  o += feed_forward(s, nn_human_SH2_SH2D3C_7, 9, 6);
  o += feed_forward(s, nn_human_SH2_SH2D3C_8, 13, 0);
  o += feed_forward(s, nn_human_SH2_SH2D3C_9, 7, 20);
  o += feed_forward(s, nn_human_SH2_SH2D3C_10, 11, 2);
  o += feed_forward(s, nn_human_SH2_SH2D3C_11, 11, 15);
  o += feed_forward(s, nn_human_SH2_SH2D3C_12, 5, 15);
  o /= 12;
  o = 0.0192057994487941+(0.066898464625402-0.0192057994487941)/(1+exp(99.9998*(0.501643-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSH2D3C\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SHC1_1, 11, 10);
  o += feed_forward(s, nn_human_SH2_SHC1_2, 9, 4);
  o += feed_forward(s, nn_human_SH2_SHC1_3, 7, 10);
  o += feed_forward(s, nn_human_SH2_SHC1_4, 11, 4);
  o += feed_forward(s, nn_human_SH2_SHC1_5, 7, 15);
  o += feed_forward(s, nn_human_SH2_SHC1_6, 13, 4);
  o += feed_forward(s, nn_human_SH2_SHC1_7, 11, 15);
  o += feed_forward(s, nn_human_SH2_SHC1_8, 7, 6);
  o += feed_forward(s, nn_human_SH2_SHC1_9, 13, 2);
  o += feed_forward(s, nn_human_SH2_SHC1_10, 9, 10);
  o += feed_forward(s, nn_human_SH2_SHC1_11, 7, 15);
  o += feed_forward(s, nn_human_SH2_SHC1_12, 11, 6);
  o /= 12;
  o = 0.0133212796771635+(0.128594182244949-0.0133212796771635)/(1+exp(100*(0.505878-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSHC1\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SHC4_1, 7, 20);
  o += feed_forward(s, nn_human_SH2_SHC4_2, 13, 15);
  o += feed_forward(s, nn_human_SH2_SHC4_3, 13, 4);
  o += feed_forward(s, nn_human_SH2_SHC4_4, 9, 15);
  o += feed_forward(s, nn_human_SH2_SHC4_5, 11, 4);
  o += feed_forward(s, nn_human_SH2_SHC4_6, 11, 15);
  o += feed_forward(s, nn_human_SH2_SHC4_7, 7, 4);
  o += feed_forward(s, nn_human_SH2_SHC4_8, 13, 6);
  o += feed_forward(s, nn_human_SH2_SHC4_9, 9, 6);
  o += feed_forward(s, nn_human_SH2_SHC4_10, 7, 2);
  o += feed_forward(s, nn_human_SH2_SHC4_11, 11, 0);
  o += feed_forward(s, nn_human_SH2_SHC4_12, 13, 15);
  o /= 12;
  o = 1.34125813862163e-07+(0.141034031413613-1.34125813862163e-07)/(1+exp(14.2293*(0.471184-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSHC4\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SHD_1, 11, 0);
  o += feed_forward(s, nn_human_SH2_SHD_2, 5, 2);
  o += feed_forward(s, nn_human_SH2_SHD_3, 5, 20);
  o += feed_forward(s, nn_human_SH2_SHD_4, 13, 20);
  o += feed_forward(s, nn_human_SH2_SHD_5, 5, 6);
  o += feed_forward(s, nn_human_SH2_SHD_6, 5, 20);
  o += feed_forward(s, nn_human_SH2_SHD_7, 7, 0);
  o += feed_forward(s, nn_human_SH2_SHD_8, 5, 2);
  o += feed_forward(s, nn_human_SH2_SHD_9, 13, 20);
  o += feed_forward(s, nn_human_SH2_SHD_10, 5, 0);
  o += feed_forward(s, nn_human_SH2_SHD_11, 5, 4);
  o += feed_forward(s, nn_human_SH2_SHD_12, 5, 20);
  o /= 12;
  o = 0.0059809942912524+(0.115781710914454-0.0059809942912524)/(1+exp(40.5056*(0.493924-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSHD\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SLAP_group_1, 13, 2);
  o += feed_forward(s, nn_human_SH2_SLAP_group_2, 13, 6);
  o += feed_forward(s, nn_human_SH2_SLAP_group_3, 13, 6);
  o += feed_forward(s, nn_human_SH2_SLAP_group_4, 7, 0);
  o += feed_forward(s, nn_human_SH2_SLAP_group_5, 11, 10);
  o += feed_forward(s, nn_human_SH2_SLAP_group_6, 5, 6);
  o += feed_forward(s, nn_human_SH2_SLAP_group_7, 7, 0);
  o += feed_forward(s, nn_human_SH2_SLAP_group_8, 13, 4);
  o += feed_forward(s, nn_human_SH2_SLAP_group_9, 13, 6);
  o += feed_forward(s, nn_human_SH2_SLAP_group_10, 7, 20);
  o += feed_forward(s, nn_human_SH2_SLAP_group_11, 7, 2);
  o += feed_forward(s, nn_human_SH2_SLAP_group_12, 13, 6);
  o /= 12;
  o = 0.00938186140991449+(0.0810704907902547-0.00938186140991449)/(1+exp(42.2825*(0.488327-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSLAP_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_1, 11, 20);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_2, 7, 4);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_3, 13, 20);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_4, 9, 20);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_5, 11, 2);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_6, 11, 15);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_7, 5, 0);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_8, 11, 10);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_9, 13, 6);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_10, 13, 15);
  o += feed_forward(s, nn_human_SH2_SOCS2_SOCS7_group_11, 11, 10);
  o /= 11;
  o = 0.0121591431154561+(0.125234711313297-0.0121591431154561)/(1+exp(28.2305*(0.489364-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSOCS2_SOCS7_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_1, 5, 20);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_2, 13, 10);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_3, 7, 6);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_4, 9, 2);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_5, 13, 6);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_6, 13, 10);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_7, 11, 6);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_8, 9, 20);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_9, 7, 10);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_10, 9, 15);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_11, 5, 20);
  o += feed_forward(s, nn_human_SH2_SOCS5_SOCS4_SOCS6_CISH_group_12, 13, 15);
  o /= 12;
  o = 0.0249422570741409+(0.251668752607426-0.0249422570741409)/(1+exp(49.3944*(0.507881-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSOCS5_SOCS4_SOCS6_CISH_group\t%.6f\t%.6f\n", o, 0.0571428571428571);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_SUPT6H_1, 9, 20);
  o += feed_forward(s, nn_human_SH2_SUPT6H_2, 7, 20);
  o += feed_forward(s, nn_human_SH2_SUPT6H_3, 7, 10);
  o += feed_forward(s, nn_human_SH2_SUPT6H_4, 9, 10);
  o += feed_forward(s, nn_human_SH2_SUPT6H_5, 7, 15);
  o += feed_forward(s, nn_human_SH2_SUPT6H_6, 7, 2);
  o += feed_forward(s, nn_human_SH2_SUPT6H_7, 7, 0);
  o += feed_forward(s, nn_human_SH2_SUPT6H_8, 7, 20);
  o += feed_forward(s, nn_human_SH2_SUPT6H_9, 9, 0);
  o += feed_forward(s, nn_human_SH2_SUPT6H_10, 11, 10);
  o += feed_forward(s, nn_human_SH2_SUPT6H_11, 11, 10);
  o += feed_forward(s, nn_human_SH2_SUPT6H_12, 9, 20);
  o /= 12;
  o = 0.0177036542511535+(0.109482235640213-0.0177036542511535)/(1+exp(62.6538*(0.519668-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSUPT6H\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Src_group_1, 11, 15);
  o += feed_forward(s, nn_human_SH2_Src_group_2, 9, 6);
  o += feed_forward(s, nn_human_SH2_Src_group_3, 7, 4);
  o += feed_forward(s, nn_human_SH2_Src_group_4, 11, 4);
  o += feed_forward(s, nn_human_SH2_Src_group_5, 11, 15);
  o += feed_forward(s, nn_human_SH2_Src_group_6, 9, 6);
  o += feed_forward(s, nn_human_SH2_Src_group_7, 13, 0);
  o += feed_forward(s, nn_human_SH2_Src_group_8, 7, 20);
  o += feed_forward(s, nn_human_SH2_Src_group_9, 13, 6);
  o += feed_forward(s, nn_human_SH2_Src_group_10, 7, 15);
  o += feed_forward(s, nn_human_SH2_Src_group_11, 5, 4);
  o += feed_forward(s, nn_human_SH2_Src_group_12, 7, 0);
  o /= 12;
  o = 0.0147906407907949+(0.348985872546597-0.0147906407907949)/(1+exp(10.4538*(0.534095-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSrc_group\t%.6f\t%.6f\n", o, 0.0808122035641769);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Srm_1, 13, 6);
  o += feed_forward(s, nn_human_SH2_Srm_2, 11, 20);
  o += feed_forward(s, nn_human_SH2_Srm_3, 13, 10);
  o += feed_forward(s, nn_human_SH2_Srm_4, 7, 4);
  o += feed_forward(s, nn_human_SH2_Srm_5, 11, 4);
  o += feed_forward(s, nn_human_SH2_Srm_6, 11, 6);
  o += feed_forward(s, nn_human_SH2_Srm_7, 5, 2);
  o += feed_forward(s, nn_human_SH2_Srm_8, 13, 15);
  o += feed_forward(s, nn_human_SH2_Srm_9, 11, 6);
  o += feed_forward(s, nn_human_SH2_Srm_10, 7, 2);
  o += feed_forward(s, nn_human_SH2_Srm_11, 7, 15);
  o += feed_forward(s, nn_human_SH2_Srm_12, 9, 15);
  o /= 12;
  o = 0.0189501131183957+(0.0782568262685913-0.0189501131183957)/(1+exp(99.9518*(0.498432-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSrm\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Syk_2_group_1, 7, 6);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_2, 7, 15);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_3, 13, 2);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_4, 9, 10);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_5, 7, 15);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_6, 7, 4);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_7, 7, 20);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_8, 9, 6);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_9, 9, 15);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_10, 9, 2);
  o += feed_forward(s, nn_human_SH2_Syk_2_group_11, 11, 4);
  o /= 11;
  o = 0.00823659923991438+(0.23029883324333-0.00823659923991438)/(1+exp(20.9046*(0.496698-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tSyk_2_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_1, 9, 10);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_2, 7, 15);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_3, 13, 20);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_4, 11, 10);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_5, 11, 15);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_6, 13, 6);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_7, 9, 0);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_8, 13, 0);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_9, 13, 20);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_10, 9, 10);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_11, 7, 2);
  o += feed_forward(s, nn_human_SH2_TENC1_TNS3_TNS1_group_12, 11, 10);
  o /= 12;
  o = 0.00289539651262628+(0.214323910399813-0.00289539651262628)/(1+exp(7.63736*(0.496404-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tTENC1_TNS3_TNS1_group\t%.6f\t%.6f\n", o, 0.0494871659305393);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_Tec_group_1, 5, 15);
  o += feed_forward(s, nn_human_SH2_Tec_group_2, 9, 0);
  o += feed_forward(s, nn_human_SH2_Tec_group_3, 7, 20);
  o += feed_forward(s, nn_human_SH2_Tec_group_4, 11, 0);
  o += feed_forward(s, nn_human_SH2_Tec_group_5, 9, 10);
  o += feed_forward(s, nn_human_SH2_Tec_group_6, 9, 15);
  o += feed_forward(s, nn_human_SH2_Tec_group_7, 9, 2);
  o += feed_forward(s, nn_human_SH2_Tec_group_8, 7, 2);
  o += feed_forward(s, nn_human_SH2_Tec_group_9, 9, 15);
  o += feed_forward(s, nn_human_SH2_Tec_group_10, 7, 20);
  o += feed_forward(s, nn_human_SH2_Tec_group_11, 7, 6);
  o += feed_forward(s, nn_human_SH2_Tec_group_12, 11, 15);
  o /= 12;
  o = 0.0181276295064537+(0.197313579156983-0.0181276295064537)/(1+exp(16.8585*(0.506623-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tTec_group\t%.6f\t%.6f\n", o, 0.063887656499994);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_1, 7, 6);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_2, 7, 0);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_3, 7, 4);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_4, 7, 10);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_5, 7, 4);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_6, 5, 4);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_7, 9, 20);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_8, 7, 20);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_9, 11, 20);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_10, 13, 20);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_11, 13, 6);
  o += feed_forward(s, nn_human_SH2_VAV1_VAV3_group_12, 7, 2);
  o /= 12;
  o = 0.00661676974609243+(0.255811318871841-0.00661676974609243)/(1+exp(78.3006*(0.507351-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tVAV1_VAV3_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_SH2_VAV2_1, 11, 4);
  o += feed_forward(s, nn_human_SH2_VAV2_2, 5, 2);
  o += feed_forward(s, nn_human_SH2_VAV2_3, 9, 10);
  o += feed_forward(s, nn_human_SH2_VAV2_4, 9, 4);
  o += feed_forward(s, nn_human_SH2_VAV2_5, 9, 15);
  o += feed_forward(s, nn_human_SH2_VAV2_6, 7, 4);
  o += feed_forward(s, nn_human_SH2_VAV2_7, 9, 0);
  o += feed_forward(s, nn_human_SH2_VAV2_8, 9, 15);
  o += feed_forward(s, nn_human_SH2_VAV2_9, 13, 10);
  o += feed_forward(s, nn_human_SH2_VAV2_10, 13, 2);
  o += feed_forward(s, nn_human_SH2_VAV2_11, 9, 6);
  o += feed_forward(s, nn_human_SH2_VAV2_12, 7, 4);
  o /= 12;
  o = 0.00340931926808764+(0.256140919961329-0.00340931926808764)/(1+exp(100*(0.504162-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tSH2\tVAV2\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTB_APPL_group_1, 9, 10);
  o += feed_forward(s, nn_human_PTB_APPL_group_2, 9, 20);
  o += feed_forward(s, nn_human_PTB_APPL_group_3, 7, 10);
  o += feed_forward(s, nn_human_PTB_APPL_group_4, 7, 2);
  o += feed_forward(s, nn_human_PTB_APPL_group_5, 13, 2);
  o += feed_forward(s, nn_human_PTB_APPL_group_6, 9, 15);
  o += feed_forward(s, nn_human_PTB_APPL_group_7, 7, 0);
  o += feed_forward(s, nn_human_PTB_APPL_group_8, 9, 15);
  o += feed_forward(s, nn_human_PTB_APPL_group_9, 5, 4);
  o += feed_forward(s, nn_human_PTB_APPL_group_10, 7, 0);
  o += feed_forward(s, nn_human_PTB_APPL_group_11, 9, 10);
  o += feed_forward(s, nn_human_PTB_APPL_group_12, 7, 2);
  o /= 12;
  o = 2.98083714771077e-10+(0.340620191086717-2.98083714771077e-10)/(1+exp(100*(0.497679-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTB\tAPPL_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTB_FRS_group_1, 13, 4);
  o += feed_forward(s, nn_human_PTB_FRS_group_2, 13, 4);
  o += feed_forward(s, nn_human_PTB_FRS_group_3, 7, 0);
  o += feed_forward(s, nn_human_PTB_FRS_group_4, 13, 4);
  o += feed_forward(s, nn_human_PTB_FRS_group_5, 7, 20);
  o += feed_forward(s, nn_human_PTB_FRS_group_6, 7, 0);
  o += feed_forward(s, nn_human_PTB_FRS_group_7, 13, 15);
  o += feed_forward(s, nn_human_PTB_FRS_group_8, 13, 20);
  o += feed_forward(s, nn_human_PTB_FRS_group_9, 7, 0);
  o += feed_forward(s, nn_human_PTB_FRS_group_10, 5, 4);
  o += feed_forward(s, nn_human_PTB_FRS_group_11, 13, 10);
  o += feed_forward(s, nn_human_PTB_FRS_group_12, 13, 0);
  o /= 12;
  o = 1.0586838967413e-06+(0.430621817443985-1.0586838967413e-06)/(1+exp(20.9758*(0.495562-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTB\tFRS_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'Y') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_PTB_SHC4_1, 13, 0);
  o += feed_forward(s, nn_human_PTB_SHC4_2, 7, 10);
  o += feed_forward(s, nn_human_PTB_SHC4_3, 7, 2);
  o += feed_forward(s, nn_human_PTB_SHC4_4, 13, 0);
  o += feed_forward(s, nn_human_PTB_SHC4_5, 7, 0);
  o += feed_forward(s, nn_human_PTB_SHC4_6, 9, 10);
  o += feed_forward(s, nn_human_PTB_SHC4_7, 9, 2);
  o += feed_forward(s, nn_human_PTB_SHC4_8, 9, 15);
  o += feed_forward(s, nn_human_PTB_SHC4_9, 9, 10);
  o += feed_forward(s, nn_human_PTB_SHC4_10, 13, 2);
  o += feed_forward(s, nn_human_PTB_SHC4_11, 7, 0);
  o += feed_forward(s, nn_human_PTB_SHC4_12, 7, 2);
  o /= 12;
  o = 3.67920081599197e-07+(0.448792845934245-3.67920081599197e-07)/(1+exp(28.03*(0.493974-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tPTB\tSHC4\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_1433_any_group_1, 13, 20);
  o += feed_forward(s, nn_human_1433_any_group_2, 13, 15);
  o += feed_forward(s, nn_human_1433_any_group_3, 7, 4);
  o += feed_forward(s, nn_human_1433_any_group_4, 9, 6);
  o += feed_forward(s, nn_human_1433_any_group_5, 13, 15);
  o += feed_forward(s, nn_human_1433_any_group_6, 7, 20);
  o += feed_forward(s, nn_human_1433_any_group_7, 13, 20);
  o += feed_forward(s, nn_human_1433_any_group_8, 11, 15);
  o += feed_forward(s, nn_human_1433_any_group_9, 13, 15);
  o += feed_forward(s, nn_human_1433_any_group_10, 5, 4);
  o += feed_forward(s, nn_human_1433_any_group_11, 7, 15);
  o += feed_forward(s, nn_human_1433_any_group_12, 11, 15);
  o /= 12;
  o = 0.0129608159538038+(0.59299965656527-0.0129608159538038)/(1+exp(29.0827*(0.531764-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\t1433\tany_group\t%.6f\t%.6f\n", o, 0.0808122035641769);
  }
}

if (c == 'S') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_WD40_any_group_1, 5, 2);
  o += feed_forward(s, nn_human_WD40_any_group_2, 5, 15);
  o += feed_forward(s, nn_human_WD40_any_group_3, 7, 2);
  o += feed_forward(s, nn_human_WD40_any_group_4, 5, 2);
  o += feed_forward(s, nn_human_WD40_any_group_5, 7, 20);
  o += feed_forward(s, nn_human_WD40_any_group_6, 5, 15);
  o += feed_forward(s, nn_human_WD40_any_group_7, 9, 10);
  o += feed_forward(s, nn_human_WD40_any_group_8, 7, 0);
  o /= 8;
  o = 0.0149729966836907+(0.091700155680349-0.0149729966836907)/(1+exp(99.9992*(0.491831-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tWD40\tany_group\t%.6f\t%.6f\n", o, 0.0404061017820884);
  }
}

if (c == 'S' || c == 'T') {
  o = 0;
  gbHas_unsupported_amino_acid = 0;
  o += feed_forward(s, nn_human_WW_any_group_1, 11, 10);
  o += feed_forward(s, nn_human_WW_any_group_2, 5, 20);
  o += feed_forward(s, nn_human_WW_any_group_3, 13, 0);
  o += feed_forward(s, nn_human_WW_any_group_4, 13, 10);
  o += feed_forward(s, nn_human_WW_any_group_5, 13, 2);
  o += feed_forward(s, nn_human_WW_any_group_6, 11, 6);
  o += feed_forward(s, nn_human_WW_any_group_7, 13, 6);
  o += feed_forward(s, nn_human_WW_any_group_8, 5, 20);
  o += feed_forward(s, nn_human_WW_any_group_9, 9, 15);
  o += feed_forward(s, nn_human_WW_any_group_10, 13, 4);
  o += feed_forward(s, nn_human_WW_any_group_11, 11, 10);
  o += feed_forward(s, nn_human_WW_any_group_12, 7, 10);
  o /= 12;
  o = 8.47863581703695e-06+(0.131893508667361-8.47863581703695e-06)/(1+exp(25.7291*(0.493791-o)));
  if (gbHas_unsupported_amino_acid == 0){
    printf("%s\t%d\t%c\t", name, *n, c);
    print_peptide(s);
    printf("\tnn\thuman\tWW\tany_group\t%.6f\t%.6f\n", o, 0.0285714285714286);
  }
}

