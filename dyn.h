#pragma once
#define N_STATES 19
void dyn(const float* x, const float* in, float* ret) {
    float subx_0 = cos(x[8]);
    float subx_1 = 0.000402*pow(subx_0, 2);
    float subx_2 = cos(x[9]);
    float subx_3 = 0.000402*pow(subx_2, 2);
    float subx_4 = sin(x[8]);
    float subx_5 = 0.000402*pow(subx_4, 2);
    float subx_6 = sin(x[9]);
    float subx_7 = 0.000402*pow(subx_6, 2);
    float subx_8 = 1.0/(subx_1 + subx_3 + subx_5 + subx_7 + 2.81223888256237);
    float subx_9 = cos(x[7]);
    float subx_10 = sin(x[7]);
    float subx_11 = subx_10*x[12] + subx_9*x[11];
    float subx_12 = -9.891512475e-5*subx_10*x[11] + 9.891512475e-5*subx_9*x[12];
    float subx_13 = -subx_11*subx_12;
    float subx_14 = -0.00241793333333333*subx_10*x[11] + 0.00241793333333333*subx_9*x[12];
    float subx_15 = -subx_11*subx_14;
    float subx_16 = -2.0e-5*subx_10*x[11] + 2.0e-5*subx_9*x[12];
    float subx_17 = -subx_11*subx_16;
    float subx_18 = x[11] + x[17];
    float subx_19 = 0.000402*subx_0*x[12] + 0.000402*subx_4*x[10];
    float subx_20 = 0.000804*x[11];
    float subx_21 = subx_20 + 0.000804*x[17];
    float subx_22 = subx_0*x[12] + subx_4*x[10];
    float subx_23 = x[11]*x[12];
    float subx_24 = -subx_18*x[12] + subx_23;
    float subx_25 = x[10]*x[11];
    float subx_26 = -subx_25;
    float subx_27 = subx_18*x[10] + subx_26;
    float subx_28 = 0.000402*subx_0*subx_24 + subx_18*subx_19 - subx_21*subx_22 - 0.000402*subx_27*subx_4;
    float subx_29 = -subx_0*subx_28;
    float subx_30 = x[11] + x[18];
    float subx_31 = 0.000402*subx_2*x[12] + 0.000402*subx_6*x[10];
    float subx_32 = subx_20 + 0.000804*x[18];
    float subx_33 = subx_2*x[12] + subx_6*x[10];
    float subx_34 = subx_23 - subx_30*x[12];
    float subx_35 = subx_26 + subx_30*x[10];
    float subx_36 = 0.000402*subx_2*subx_34 + subx_30*subx_31 - subx_32*subx_33 - 0.000402*subx_35*subx_6;
    float subx_37 = -subx_2*subx_36;
    float subx_38 = subx_0*x[10] - subx_4*x[12];
    float subx_39 = 0.000402*subx_0*x[10] - 0.000402*subx_4*x[12];
    float subx_40 = 0.000402*subx_0*subx_27 - subx_18*subx_39 + subx_21*subx_38 + 0.000402*subx_24*subx_4;
    float subx_41 = -subx_4*subx_40;
    float subx_42 = subx_2*x[10] - subx_6*x[12];
    float subx_43 = 0.000402*subx_2*x[10] - 0.000402*subx_6*x[12];
    float subx_44 = 0.000402*subx_2*subx_35 - subx_30*subx_43 + subx_32*subx_42 + 0.000402*subx_34*subx_6;
    float subx_45 = -subx_44*subx_6;
    float subx_46 = 0.0904682787687477*subx_10*subx_8 - 0.03225178*subx_10;
    float subx_47 = pow(subx_9, 2);
    float subx_48 = pow(subx_10, 2);
    float subx_49 = 2.80311499304904*subx_47 - 0.0010401773131684*subx_48*subx_8 + 0.00298411727808333*subx_48 + 0.002888;
    float subx_50 = 1.0/subx_49;
    float subx_51 = 0.0010401773131684*subx_10*subx_8*subx_9 + 2.80013087577096*subx_10*subx_9;
    float subx_52 = -subx_46*subx_50*subx_51 - 0.0904682787687477*subx_8*subx_9 + 0.03225178*subx_9;
    float subx_53 = pow(subx_51, 2);
    float subx_54 = 1.0/(subx_1 + subx_3 - 0.0010401773131684*subx_47*subx_8 + 0.00298411727808333*subx_47 + 2.80311499304904*subx_48 + subx_5 - subx_50*subx_53 + subx_7 + 0.00717575833333333);
    float subx_55 = 2*x[0]*x[3];
    float subx_56 = 2*x[1]*x[2];
    float subx_57 = -subx_55 + subx_56;
    float subx_58 = 1.82685*subx_57*subx_9;
    float subx_59 = 2*x[0]*x[2];
    float subx_60 = 2*x[1]*x[3];
    float subx_61 = subx_59 + subx_60;
    float subx_62 = 1.82685*subx_10*subx_61;
    float subx_63 = subx_58 + subx_62;
    float subx_64 = subx_10*subx_61 + subx_57*subx_9;
    float subx_65 = -subx_10*subx_57 + subx_61*subx_9;
    float subx_66 = pow(x[0], 2);
    float subx_67 = pow(x[1], 2);
    float subx_68 = pow(x[2], 2);
    float subx_69 = -subx_68;
    float subx_70 = pow(x[3], 2);
    float subx_71 = -subx_70;
    float subx_72 = subx_66 + subx_67 + subx_69 + subx_71;
    float subx_73 = 0.03225178*subx_10*subx_63*subx_8 - 0.018916*subx_10*subx_64 - 0.018916*subx_65*subx_9 - 1.82685*subx_72*subx_9;
    float subx_74 = -0.018916*subx_10*subx_65 - 1.82685*subx_10*subx_72 - subx_50*subx_51*subx_73 - 0.03225178*subx_63*subx_8*subx_9 + 0.018916*subx_64*subx_9;
    float subx_75 = 1.0/(-subx_50*pow(subx_73, 2) - subx_54*pow(subx_74, 2) - pow(subx_63, 2)*subx_8 + 2.39);
    float subx_76 = -subx_46*subx_50*subx_73 - subx_52*subx_54*subx_74 + subx_58 + subx_62 - 2.80506312422904*subx_63*subx_8;
    float subx_77 = 2*x[0]*x[1];
    float subx_78 = -subx_77;
    float subx_79 = 2*x[2]*x[3];
    float subx_80 = subx_78 + subx_79;
    float subx_81 = 1.82685*subx_10*subx_80;
    float subx_82 = -subx_67;
    float subx_83 = subx_66 + subx_68 + subx_71 + subx_82;
    float subx_84 = 1.82685*subx_83*subx_9;
    float subx_85 = subx_81 + subx_84;
    float subx_86 = -subx_10*subx_83 + subx_80*subx_9;
    float subx_87 = subx_10*subx_80 + subx_83*subx_9;
    float subx_88 = subx_55 + subx_56;
    float subx_89 = 0.03225178*subx_10*subx_8*subx_85 - 0.018916*subx_10*subx_87 - 0.018916*subx_86*subx_9 - 1.82685*subx_88*subx_9;
    float subx_90 = -0.018916*subx_10*subx_86 - 1.82685*subx_10*subx_88 - subx_50*subx_51*subx_89 - 0.03225178*subx_8*subx_85*subx_9 + 0.018916*subx_87*subx_9;
    float subx_91 = -subx_50*subx_73*subx_89 - subx_54*subx_74*subx_90 - subx_63*subx_8*subx_85;
    float subx_92 = 1.0/(-subx_50*pow(subx_89, 2) - subx_54*pow(subx_90, 2) - subx_75*pow(subx_91, 2) - subx_8*pow(subx_85, 2) + 2.39);
    float subx_93 = -subx_46*subx_50*subx_89 - subx_52*subx_54*subx_90 - subx_75*subx_76*subx_91 - 2.80506312422904*subx_8*subx_85 + subx_81 + subx_84;
    float subx_94 = subx_77 + subx_79;
    float subx_95 = 1.82685*subx_9*subx_94;
    float subx_96 = subx_66 + subx_69 + subx_70 + subx_82;
    float subx_97 = 1.82685*subx_10*subx_96;
    float subx_98 = subx_95 + subx_97;
    float subx_99 = subx_10*subx_96 + subx_9*subx_94;
    float subx_100 = -subx_10*subx_94 + subx_9*subx_96;
    float subx_101 = -subx_59 + subx_60;
    float subx_102 = 0.03225178*subx_10*subx_8*subx_98 - 0.018916*subx_10*subx_99 - 0.018916*subx_100*subx_9 - 1.82685*subx_101*subx_9;
    float subx_103 = -0.018916*subx_10*subx_100 - 1.82685*subx_10*subx_101 - subx_102*subx_50*subx_51 - 0.03225178*subx_8*subx_9*subx_98 + 0.018916*subx_9*subx_99;
    float subx_104 = -subx_102*subx_50*subx_73 - subx_103*subx_54*subx_74 - subx_63*subx_8*subx_98;
    float subx_105 = -subx_102*subx_50*subx_89 - subx_103*subx_54*subx_90 - subx_104*subx_75*subx_91 - subx_8*subx_85*subx_98;
    float subx_106 = 1.0/(-pow(subx_102, 2)*subx_50 - pow(subx_103, 2)*subx_54 - pow(subx_104, 2)*subx_75 - pow(subx_105, 2)*subx_92 - subx_8*pow(subx_98, 2) + 2.39);
    float subx_107 = -subx_102*subx_46*subx_50 - subx_103*subx_52*subx_54 - subx_104*subx_75*subx_76 - subx_105*subx_92*subx_93 - 2.80506312422904*subx_8*subx_98 + subx_95 + subx_97;
    float subx_108 = 1.0/(-subx_106*pow(subx_107, 2) - pow(subx_46, 2)*subx_50 - pow(subx_52, 2)*subx_54 - subx_75*pow(subx_76, 2) - 7.86837913090959*subx_8 - subx_92*pow(subx_93, 2) + 2.80506312422904);
    float subx_109 = -subx_22*subx_39;
    float subx_110 = -subx_33*subx_43;
    float subx_111 = 0.00937333333333334*x[10] + 0.00937333333333334*x[16];
    float subx_112 = -subx_10*x[11] + subx_9*x[12];
    float subx_113 = x[10] + x[16];
    float subx_114 = -subx_113*x[11] + subx_25;
    float subx_115 = subx_113*x[12] - x[10]*x[12];
    float subx_116 = 0.00937333333333334*subx_10*subx_114 + subx_111*subx_112 - subx_113*subx_16 + 0.00937333333333334*subx_115*subx_9;
    float subx_117 = -subx_116*subx_9;
    float subx_118 = 0.00657333333333333*x[10] + 0.00657333333333333*x[16];
    float subx_119 = 0.00417793333333333*subx_10*subx_114 + subx_112*subx_118 - subx_113*subx_14 + 0.00417793333333333*subx_115*subx_9;
    float subx_120 = -subx_119*subx_9;
    float subx_121 = 0.108888707562375*x[10] + 0.108888707562375*x[16];
    float subx_122 = 0.108888707562375*subx_10*subx_114 + subx_112*subx_121 - subx_113*subx_12 + 0.108888707562375*subx_115*subx_9;
    float subx_123 = -subx_122*subx_9;
    float subx_124 = subx_19*subx_38;
    float subx_125 = subx_31*subx_42;
    float subx_126 = 0.00417793333333333*subx_10*x[12] + 0.00417793333333333*subx_9*x[11];
    float subx_127 = -0.00241793333333333*subx_10*subx_115 - subx_11*subx_118 + subx_113*subx_126 + 0.00241793333333333*subx_114*subx_9;
    float subx_128 = 0.00937333333333334*subx_10*x[12] + 0.00937333333333334*subx_9*x[11];
    float subx_129 = -2.0e-5*subx_10*subx_115 - subx_11*subx_111 + subx_113*subx_128 + 2.0e-5*subx_114*subx_9;
    float subx_130 = 0.108888707562375*subx_10*x[12] + 0.108888707562375*subx_9*x[11];
    float subx_131 = -9.891512475e-5*subx_10*subx_115 - subx_11*subx_121 + subx_113*subx_130 + 9.891512475e-5*subx_114*subx_9;
    float subx_132 = -0.023645*subx_10*x[11] + 0.023645*subx_9*x[12] + 1.705*x[10] + 1.705*x[16];
    float subx_133 = -1.705*subx_10*x[12] - 1.705*subx_9*x[11];
    float subx_134 = 0.018916*subx_10*x[11]*x[16] - 0.8*subx_11*subx_133 + 0.8*subx_113*subx_132 - 0.018916*subx_9*x[12]*x[16];
    float subx_135 = -0.023645*subx_10*x[12] - 0.023645*subx_9*x[11];
    float subx_136 = -0.018916*subx_10*x[12]*x[16] + 0.8*subx_112*subx_133 - 0.8*subx_113*subx_135 - 0.018916*subx_9*x[11]*x[16];
    float subx_137 = 1.364*subx_10*x[11]*x[16] + 0.8*subx_11*subx_135 - 0.8*subx_112*subx_132 - 1.364*subx_9*x[12]*x[16];
    float subx_138 = 0.915*x[10] + 0.915*x[16];
    float subx_139 = 0.35685*subx_10*x[11]*x[16] - 0.39*subx_112*subx_138 - 0.35685*subx_9*x[12]*x[16];
    float subx_140 = 0.265*x[10] + 0.265*x[16];
    float subx_141 = 0.106*subx_10*x[11]*x[16] - 0.4*subx_112*subx_140 - 0.106*subx_9*x[12]*x[16];
    float subx_142 = 7.84524*subx_10*subx_96 + 7.84524*subx_9*subx_94;
    float subx_143 = -7.84524*subx_10*subx_94 + 7.84524*subx_9*subx_96;
    float subx_144 = -15.69048*x[0]*x[2] + 15.69048*x[1]*x[3];
    float subx_145 = -7.649109*x[0]*x[2] + 7.649109*x[1]*x[3];
    float subx_146 = -7.84524*x[0]*x[2] + 7.84524*x[1]*x[3];
    float subx_147 = subx_78 - subx_79;
    float subx_148 = pow(pow(subx_147, 2) + 2*subx_147*subx_94 + 1, -1.0L/2.0L);
    float subx_149 = 0.092*subx_147*subx_148;
    float subx_150 = subx_149 + 0.11375;
    float subx_151 = 920.0*subx_148;
    float subx_152 = 10000.0*x[6];
    float subx_153 = 0.5*erf(10000.0*subx_150*subx_94 + subx_151 + subx_152) + 0.5;
    float subx_154 = pow(M_PI, 2);
    float subx_155 = 0.092*subx_148;
    float subx_156 = 0.11375*x[10];
    float subx_157 = 0.092*subx_147*subx_148*x[10];
    float subx_158 = subx_156 + subx_157;
    float subx_159 = 0.11375*x[12];
    float subx_160 = -0.092*subx_147*subx_148*x[12];
    float subx_161 = -subx_159 + subx_160;
    float subx_162 = -4302.0*subx_154*(subx_150*subx_94 + subx_155 + x[6]) - 43.02*M_PI*(subx_101*subx_161 + subx_158*subx_96 + x[15]);
    float subx_163 = subx_153*subx_162;
    float subx_164 = -subx_163;
    float subx_165 = ((subx_164 < 0) ? (
   0
)
: (
   subx_164
));
    float subx_166 = subx_88*x[10];
    float subx_167 = subx_80*x[12];
    float subx_168 = 0.092*subx_148*(subx_166 + subx_167 + subx_30*subx_83) + subx_158*subx_61 + subx_161*subx_72 + x[13];
    float subx_169 = subx_72*x[10];
    float subx_170 = subx_61*x[12];
    float subx_171 = -0.092*subx_148*(subx_169 + subx_170 + subx_30*subx_57) + subx_158*subx_80 + subx_161*subx_88 + x[14];
    float subx_172 = sqrt(pow(fabs(subx_168), 2) + pow(fabs(subx_171), 2));
    float subx_173 = 1.0/subx_172;
    float subx_174 = ((subx_172 > 0.0) ? (
   subx_171*subx_173
)
: (
   0.0
));
    float subx_175 = ((subx_172 < 0.0001) ? (
   12000.0*subx_172
)
: ((subx_172 >= 0.0001 && subx_172 < 0.0002) ? (
   -2000.0*subx_172 + 1.4
)
: (
   1.0
)));
    float subx_176 = ((subx_172 > 0.0) ? (
   subx_168*subx_173
)
: (
   1.0
));
    float subx_177 = subx_149 - 0.11375;
    float subx_178 = 0.5*erf(subx_151 + subx_152 + 10000.0*subx_177*subx_94) + 0.5;
    float subx_179 = subx_159 + subx_160;
    float subx_180 = -subx_156 + subx_157;
    float subx_181 = -4302.0*subx_154*(subx_155 + subx_177*subx_94 + x[6]) - 43.02*M_PI*(subx_101*subx_179 + subx_180*subx_96 + x[15]);
    float subx_182 = subx_178*subx_181;
    float subx_183 = -subx_182;
    float subx_184 = ((subx_183 < 0) ? (
   0
)
: (
   subx_183
));
    float subx_185 = 0.092*subx_148*(subx_166 + subx_167 + subx_18*subx_83) + subx_179*subx_72 + subx_180*subx_61 + x[13];
    float subx_186 = -0.092*subx_148*(subx_169 + subx_170 + subx_18*subx_57) + subx_179*subx_88 + subx_180*subx_80 + x[14];
    float subx_187 = sqrt(pow(fabs(subx_185), 2) + pow(fabs(subx_186), 2));
    float subx_188 = 1.0/subx_187;
    float subx_189 = ((subx_187 > 0.0) ? (
   subx_186*subx_188
)
: (
   0.0
));
    float subx_190 = ((subx_187 < 0.0001) ? (
   12000.0*subx_187
)
: ((subx_187 >= 0.0001 && subx_187 < 0.0002) ? (
   -2000.0*subx_187 + 1.4
)
: (
   1.0
)));
    float subx_191 = ((subx_187 > 0.0) ? (
   subx_185*subx_188
)
: (
   1.0
));
    float subx_192 = subx_112*subx_126;
    float subx_193 = subx_112*subx_128;
    float subx_194 = subx_112*subx_130;
    float subx_195 = 17.9150958675*subx_9*subx_94;
    float subx_196 = 17.9150958675*subx_10*subx_96;
    float subx_197 = 1.364*subx_113*subx_135;
    float subx_198 = -0.915*subx_10*x[12] - 0.915*subx_9*x[11];
    float subx_199 = -0.35685*subx_112*subx_198;
    float subx_200 = -0.265*subx_10*x[12] - 0.265*subx_9*x[11];
    float subx_201 = -0.106*subx_112*subx_200;
    float subx_202 = -1.364*subx_112*subx_133;
    float subx_203 = 0.03225178*subx_9*x[11]*x[16];
    float subx_204 = 0.03225178*subx_10*x[12]*x[16];
    float subx_205 = 0.5*erf(18300.0*subx_10*subx_94 + subx_152 - 18300.0*subx_9*subx_96) + 0.5;
    float subx_206 = -4302.0*subx_154*(1.83*subx_10*subx_94 - 1.83*subx_9*subx_96 + x[6]) - 43.02*M_PI*(subx_101*(-1.83*subx_10*x[12] - 1.83*subx_9*x[11]) + subx_99*(1.83*x[10] + 1.83*x[16]) + x[15]);
    float subx_207 = 1.83*subx_205*subx_206*subx_99;
    float subx_208 = 0.092*subx_148*subx_165*subx_174*subx_175;
    float subx_209 = 0.092*subx_148*subx_184*subx_189*subx_190;
    float subx_210 = 0.092*subx_148*subx_165*subx_175;
    float subx_211 = 0.092*subx_148*subx_184*subx_190;
    float subx_212 = subx_150*(subx_153*subx_162*subx_96 - subx_165*subx_174*subx_175*subx_80 - subx_165*subx_175*subx_176*subx_61) - subx_176*subx_210*subx_88 + subx_177*(subx_178*subx_181*subx_96 - subx_184*subx_189*subx_190*subx_80 - subx_184*subx_190*subx_191*subx_61) - subx_191*subx_211*subx_88 + subx_192 + subx_193 + subx_194 + subx_195 + subx_196 + subx_197 + subx_199 + subx_201 + subx_202 + subx_203 + subx_204 + subx_207 + subx_208*subx_72 + subx_209*subx_72 - 0.00589575833333333*x[11]*x[12];
    float subx_213 = subx_13 + subx_15 + subx_17 + subx_212 + subx_29 + subx_37 + subx_41 + subx_45;
    float subx_214 = 1.83*subx_101*subx_205*subx_206;
    float subx_215 = subx_208*subx_57;
    float subx_216 = subx_189*subx_211*subx_57;
    float subx_217 = 0.092*subx_148*subx_165*subx_175*subx_176;
    float subx_218 = -subx_217*subx_83;
    float subx_219 = 0.092*subx_148*subx_184*subx_190*subx_191;
    float subx_220 = -subx_219*subx_83;
    float subx_221 = subx_10*subx_127 + subx_10*subx_129 + subx_10*subx_131 + 0.023645*subx_10*subx_136 - 0.023645*subx_10*subx_142 + 0.03225178*subx_10*subx_213*subx_8 + subx_124 + subx_125 + 0.023645*subx_134*subx_9 + 1.705*subx_137*subx_9 + 0.915*subx_139*subx_9 + 0.265*subx_141*subx_9 - 0.023645*subx_143*subx_9 - 1.705*subx_144*subx_9 - 0.915*subx_145*subx_9 - 0.265*subx_146*subx_9 - subx_214*subx_9 + subx_215 + subx_216 + subx_218 + subx_220;
    float subx_222 = subx_109 + subx_110 + subx_117 + subx_120 + subx_123 + subx_221;
    float subx_223 = -subx_222*subx_46*subx_50;
    float subx_224 = -0.39*subx_112*subx_198*subx_64 - 0.4*subx_112*subx_200*subx_64;
    float subx_225 = 0.02275*pow(x[10], 2);
    float subx_226 = 0.02275*pow(x[12], 2);
    float subx_227 = subx_225 + subx_226;
    float subx_228 = -subx_227*subx_57;
    float subx_229 = -subx_225 - subx_226;
    float subx_230 = -subx_229*subx_57;
    float subx_231 = -subx_136*subx_64;
    float subx_232 = -0.39*subx_11*subx_198 + 0.39*subx_113*subx_138;
    float subx_233 = -subx_232*subx_65;
    float subx_234 = -0.4*subx_11*subx_200 + 0.4*subx_113*subx_140;
    float subx_235 = -subx_234*subx_65;
    float subx_236 = -subx_134*subx_65;
    float subx_237 = -subx_139*subx_72;
    float subx_238 = -subx_141*subx_72;
    float subx_239 = -subx_137*subx_72;
    float subx_240 = -subx_222*subx_50*subx_73;
    float subx_241 = -subx_213*subx_63*subx_8;
    float subx_242 = -subx_127*subx_9;
    float subx_243 = -subx_129*subx_9;
    float subx_244 = -subx_131*subx_9;
    float subx_245 = -subx_0*subx_40;
    float subx_246 = -subx_2*subx_44;
    float subx_247 = -subx_10*subx_116;
    float subx_248 = -subx_10*subx_119;
    float subx_249 = -subx_10*subx_122;
    float subx_250 = -subx_222*subx_50*subx_51;
    float subx_251 = -subx_149;
    float subx_252 = 0.023645*subx_10*subx_134 + 1.705*subx_10*subx_137 + 0.915*subx_10*subx_139 + 0.265*subx_10*subx_141 - 0.023645*subx_10*subx_143 - 1.705*subx_10*subx_144 - 0.915*subx_10*subx_145 - 0.265*subx_10*subx_146 - subx_10*subx_214 - 0.023645*subx_136*subx_9 + 0.023645*subx_142*subx_9 + subx_174*subx_210*subx_61 + subx_209*subx_61 - 0.03225178*subx_213*subx_8*subx_9 - subx_217*subx_80 - subx_219*subx_80 + subx_28*subx_4 + subx_36*subx_6 + 0.00589575833333333*x[10]*x[11] + (subx_251 - 0.11375)*(subx_101*subx_153*subx_162 - subx_165*subx_174*subx_175*subx_88 - subx_165*subx_175*subx_176*subx_72) + (subx_251 + 0.11375)*(subx_101*subx_178*subx_181 - subx_184*subx_189*subx_190*subx_88 - subx_184*subx_190*subx_191*subx_72);
    float subx_253 = subx_242 + subx_243 + subx_244 + subx_245 + subx_246 + subx_247 + subx_248 + subx_249 + subx_250 + subx_252;
    float subx_254 = -subx_253*subx_54*subx_74;
    float subx_255 = -subx_165*subx_175*subx_176;
    float subx_256 = -subx_184*subx_190*subx_191;
    float subx_257 = subx_224 + subx_228 + subx_230 + subx_231 + subx_233 + subx_235 + subx_236 + subx_237 + subx_238 + subx_239 + subx_240 + subx_241 + subx_254 + subx_255 + subx_256;
    float subx_258 = -subx_257*subx_75*subx_76;
    float subx_259 = -0.39*subx_112*subx_198*subx_87 - 0.4*subx_112*subx_200*subx_87;
    float subx_260 = -subx_227*subx_83;
    float subx_261 = -subx_229*subx_83;
    float subx_262 = -subx_232*subx_86;
    float subx_263 = -subx_234*subx_86;
    float subx_264 = -subx_134*subx_86;
    float subx_265 = -subx_136*subx_87;
    float subx_266 = -subx_139*subx_88;
    float subx_267 = -subx_141*subx_88;
    float subx_268 = -subx_137*subx_88;
    float subx_269 = -subx_222*subx_50*subx_89;
    float subx_270 = -subx_257*subx_75*subx_91;
    float subx_271 = -subx_213*subx_8*subx_85;
    float subx_272 = -subx_253*subx_54*subx_90;
    float subx_273 = -subx_165*subx_174*subx_175;
    float subx_274 = -subx_184*subx_189*subx_190;
    float subx_275 = subx_259 + subx_260 + subx_261 + subx_262 + subx_263 + subx_264 + subx_265 + subx_266 + subx_267 + subx_268 + subx_269 + subx_270 + subx_271 + subx_272 + subx_273 + subx_274;
    float subx_276 = -subx_275*subx_92*subx_93;
    float subx_277 = -0.39*subx_112*subx_198*subx_99 - 0.4*subx_112*subx_200*subx_99 + subx_163 + subx_182 + subx_205*subx_206 + 23.4376545;
    float subx_278 = -subx_227*subx_94;
    float subx_279 = -subx_229*subx_94;
    float subx_280 = -subx_136*subx_99;
    float subx_281 = -subx_100*subx_232;
    float subx_282 = -subx_100*subx_234;
    float subx_283 = -subx_100*subx_134;
    float subx_284 = -subx_101*subx_139;
    float subx_285 = -subx_101*subx_141;
    float subx_286 = -subx_101*subx_137;
    float subx_287 = -subx_102*subx_222*subx_50;
    float subx_288 = -subx_104*subx_257*subx_75;
    float subx_289 = -subx_213*subx_8*subx_98;
    float subx_290 = -subx_105*subx_275*subx_92;
    float subx_291 = -subx_103*subx_253*subx_54;
    float subx_292 = subx_277 + subx_278 + subx_279 + subx_280 + subx_281 + subx_282 + subx_283 + subx_284 + subx_285 + subx_286 + subx_287 + subx_288 + subx_289 + subx_290 + subx_291;
    float subx_293 = -subx_106*subx_107*subx_292;
    float subx_294 = -6.46416e-7*subx_50;
    float subx_295 = 0.000804*subx_50*subx_51*subx_54*subx_74 - 0.000804*subx_50*subx_73;
    float subx_296 = -pow(subx_295, 2)*subx_75;
    float subx_297 = -subx_295*subx_75*subx_91 + 0.000804*subx_50*subx_51*subx_54*subx_90 - 0.000804*subx_50*subx_89;
    float subx_298 = -pow(subx_297, 2)*subx_92;
    float subx_299 = -0.000804*subx_102*subx_50 + 0.000804*subx_103*subx_50*subx_51*subx_54 - subx_104*subx_295*subx_75 - subx_105*subx_297*subx_92;
    float subx_300 = -subx_106*pow(subx_299, 2);
    float subx_301 = -subx_106*subx_107*subx_299 - subx_295*subx_75*subx_76 - subx_297*subx_92*subx_93 - 0.000804*subx_46*subx_50 + 0.000804*subx_50*subx_51*subx_52*subx_54;
    float subx_302 = -subx_108*pow(subx_301, 2);
    float subx_303 = -6.46416e-7*subx_53*subx_54/pow(subx_49, 2);
    float subx_304 = 1.0/(subx_294 + subx_296 + subx_298 + subx_300 + subx_302 + subx_303 + 0.000804);
    float subx_305 = -0.000804*subx_222*subx_50;
    float subx_306 = -subx_257*subx_295*subx_75;
    float subx_307 = -subx_275*subx_297*subx_92;
    float subx_308 = -subx_106*subx_292*subx_299;
    float subx_309 = -subx_253*subx_52*subx_54;
    float subx_310 = -25.2455681180614*subx_154*x[7] + subx_192 + subx_193 + subx_194 + subx_195 + subx_196 + subx_197 + subx_199 + subx_201 + subx_202 + subx_203 + subx_204 + subx_207 - 2.80506312422904*subx_213*subx_8 - 8.41518937268713*M_PI*x[16];
    float subx_311 = -subx_108*subx_301*(subx_13 + subx_15 + subx_17 + subx_223 + subx_258 + subx_276 + subx_293 + subx_309 + subx_310);
    float subx_312 = subx_294 + subx_296 + subx_298 + subx_300 + subx_302 + subx_303;
    float subx_313 = 1.0/(subx_294 + subx_296 + subx_298 + subx_300 + subx_302 + subx_303 - subx_304*pow(subx_312, 2) + 0.000804);
    float subx_314 = 0.000804*subx_253*subx_50*subx_51*subx_54;
    float subx_315 = in[0] + subx_110 + subx_125 + subx_215 + subx_218 - subx_304*subx_312*(in[1] + subx_109 + subx_124 + subx_216 + subx_220 + subx_305 + subx_306 + subx_307 + subx_308 + subx_311 + subx_314) + subx_305 + subx_306 + subx_307 + subx_308 + subx_311 + subx_314;
    float subx_316 = in[1] + subx_109 + subx_124 + subx_216 + subx_220 + subx_305 + subx_306 + subx_307 + subx_308 + subx_311 - subx_312*subx_313*subx_315 + subx_314;
    float subx_317 = subx_13 + subx_15 + subx_17 + subx_223 + subx_258 + subx_276 + subx_293 - subx_301*subx_304*subx_316 - subx_301*subx_313*subx_315 + subx_309 + subx_310;
    float subx_318 = -subx_107*subx_108*subx_317 + subx_277 + subx_278 + subx_279 + subx_280 + subx_281 + subx_282 + subx_283 + subx_284 + subx_285 + subx_286 + subx_287 + subx_288 + subx_289 + subx_290 + subx_291 - subx_299*subx_304*subx_316 - subx_299*subx_313*subx_315;
    float subx_319 = -subx_105*subx_106*subx_318 - subx_108*subx_317*subx_93 + subx_259 + subx_260 + subx_261 + subx_262 + subx_263 + subx_264 + subx_265 + subx_266 + subx_267 + subx_268 + subx_269 + subx_270 + subx_271 + subx_272 + subx_273 + subx_274 - subx_297*subx_304*subx_316 - subx_297*subx_313*subx_315;
    float subx_320 = -subx_104*subx_106*subx_318 - subx_108*subx_317*subx_76 + subx_224 + subx_228 + subx_230 + subx_231 + subx_233 + subx_235 + subx_236 + subx_237 + subx_238 + subx_239 + subx_240 + subx_241 + subx_254 + subx_255 + subx_256 - subx_295*subx_304*subx_316 - subx_295*subx_313*subx_315 - subx_319*subx_91*subx_92;
    float subx_321 = -subx_103*subx_106*subx_318 - subx_108*subx_317*subx_52 + subx_242 + subx_243 + subx_244 + subx_245 + subx_246 + subx_247 + subx_248 + subx_249 + subx_250 + subx_252 + 0.000804*subx_304*subx_316*subx_50*subx_51 + 0.000804*subx_313*subx_315*subx_50*subx_51 - subx_319*subx_90*subx_92 - subx_320*subx_74*subx_75;
    float subx_322 = -subx_102*subx_106*subx_318 - subx_108*subx_317*subx_46 + subx_109 + subx_110 + subx_117 + subx_120 + subx_123 + subx_221 - 0.000804*subx_304*subx_316 - 0.000804*subx_313*subx_315 - subx_319*subx_89*subx_92 - subx_320*subx_73*subx_75 - subx_321*subx_51*subx_54;

ret[0] = -0.5*x[10]*x[1] - 0.5*x[11]*x[2] - 0.5*x[12]*x[3]
ret[1] = 0.5*x[0]*x[10] - 0.5*x[11]*x[3] + 0.5*x[12]*x[2]
ret[2] = 0.5*x[0]*x[11] + 0.5*x[10]*x[3] - 0.5*x[12]*x[1]
ret[3] = 0.5*x[0]*x[12] - 0.5*x[10]*x[2] + 0.5*x[11]*x[1]
ret[4] = x[13]
ret[5] = x[14]
ret[6] = x[15]
ret[7] = x[16]
ret[8] = x[17]
ret[9] = x[18]
ret[10] = subx_8*(0.03225178*subx_10*subx_322*subx_50 - subx_106*subx_318*subx_98 - 2.80506312422904*subx_108*subx_317 + subx_13 + subx_15 + subx_17 + subx_212 + subx_29 - subx_319*subx_85*subx_92 - subx_320*subx_63*subx_75 - 0.03225178*subx_321*subx_54*subx_9 + subx_37 + subx_41 + subx_45)
ret[11] = subx_322*subx_50
ret[12] = subx_321*subx_54
ret[13] = subx_320*subx_75
ret[14] = subx_319*subx_92
ret[15] = subx_106*subx_318
ret[16] = subx_108*subx_317
ret[17] = subx_304*subx_316
ret[18] = subx_313*subx_315
}
