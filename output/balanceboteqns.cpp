#include "balanceboteqns.h"

Eigen::Matrix<double,19,19> get_mm(Eigen::Matrix<double,19,1> states, Eigen::Matrix<double,2,1> inputs) {
    double subx_0 = states(8,0);
    double subx_1 = cos(subx_0);
    double subx_2 = states(9,0);
    double subx_3 = cos(subx_2);
    double subx_4 = sin(subx_0);
    double subx_5 = sin(subx_2);
    double subx_6 = 0.000402*pow(subx_1, 2) + 0.000402*pow(subx_3, 2) + 0.000402*pow(subx_4, 2) + 0.000402*pow(subx_5, 2);
    double subx_7 = states(7,0);
    double subx_8 = sin(subx_7);
    double subx_9 = -0.03225178*subx_8;
    double subx_10 = cos(subx_7);
    double subx_11 = 0.03225178*subx_10;
    double subx_12 = states(3,0);
    double subx_13 = states(0,0);
    double subx_14 = 2*subx_13;
    double subx_15 = subx_12*subx_14;
    double subx_16 = states(1,0);
    double subx_17 = states(2,0);
    double subx_18 = 2*subx_17;
    double subx_19 = subx_16*subx_18;
    double subx_20 = -subx_15 + subx_19;
    double subx_21 = subx_10*subx_20;
    double subx_22 = subx_13*subx_18;
    double subx_23 = subx_12*subx_16;
    double subx_24 = 2*subx_23;
    double subx_25 = subx_22 + subx_24;
    double subx_26 = subx_25*subx_8;
    double subx_27 = 1.82685*subx_21 + 1.82685*subx_26;
    double subx_28 = subx_14*subx_16;
    double subx_29 = -subx_28;
    double subx_30 = subx_12*subx_18;
    double subx_31 = subx_29 + subx_30;
    double subx_32 = subx_31*subx_8;
    double subx_33 = pow(subx_17, 2);
    double subx_34 = pow(subx_12, 2);
    double subx_35 = -subx_34;
    double subx_36 = pow(subx_13, 2);
    double subx_37 = pow(subx_16, 2);
    double subx_38 = subx_36 - subx_37;
    double subx_39 = subx_33 + subx_35 + subx_38;
    double subx_40 = subx_10*subx_39;
    double subx_41 = 1.82685*subx_32 + 1.82685*subx_40;
    double subx_42 = subx_28 + subx_30;
    double subx_43 = subx_10*subx_42;
    double subx_44 = -subx_33;
    double subx_45 = subx_34 + subx_38 + subx_44;
    double subx_46 = subx_45*subx_8;
    double subx_47 = 1.82685*subx_43 + 1.82685*subx_46;
    double subx_48 = pow(subx_10, 2);
    double subx_49 = pow(subx_8, 2);
    double subx_50 = 2.80013087577096*subx_10*subx_8;
    double subx_51 = subx_21 + subx_26;
    double subx_52 = 0.018916*subx_8;
    double subx_53 = subx_10*subx_25 - subx_20*subx_8;
    double subx_54 = 0.018916*subx_10;
    double subx_55 = -1.82685*subx_33 - 1.82685*subx_34 + 1.82685*subx_36 + 1.82685*subx_37;
    double subx_56 = -subx_10*subx_55 - subx_51*subx_52 - subx_53*subx_54;
    double subx_57 = subx_10*subx_31 - subx_39*subx_8;
    double subx_58 = subx_32 + subx_40;
    double subx_59 = subx_15 + subx_19;
    double subx_60 = 1.82685*subx_10;
    double subx_61 = -subx_52*subx_58 - subx_54*subx_57 - subx_59*subx_60;
    double subx_62 = subx_43 + subx_46;
    double subx_63 = subx_10*subx_45;
    double subx_64 = subx_42*subx_8;
    double subx_65 = subx_63 - subx_64;
    double subx_66 = -subx_22 + subx_24;
    double subx_67 = -subx_52*subx_62 - subx_54*subx_65 - subx_60*subx_66;
    double subx_68 = subx_51*subx_54 - subx_52*subx_53 - subx_55*subx_8;
    double subx_69 = 1.82685*subx_8;
    double subx_70 = -subx_52*subx_57 + subx_54*subx_58 - subx_59*subx_69;
    double subx_71 = -subx_52*subx_65 + subx_54*subx_62 - subx_66*subx_69;
    double subx_72 = states(10,0);
    double subx_73 = 0.5*subx_72;
    double subx_74 = states(11,0);
    double subx_75 = 0.5*subx_17;
    double subx_76 = states(12,0);
    double subx_77 = 0.5*subx_76;
    double subx_78 = 0.5*subx_74;
    double subx_79 = states(13,0);
    double subx_80 = states(14,0);
    double subx_81 = states(15,0);
    double subx_82 = states(16,0);
    double subx_83 = states(17,0);
    double subx_84 = states(18,0);
    double subx_85 = subx_29 - subx_30;
    double subx_86 = pow(2*subx_42*subx_85 + pow(subx_85, 2) + 1, -1.0L/2.0L);
    double subx_87 = 0.092*subx_86;
    double subx_88 = subx_85*subx_87;
    double subx_89 = subx_88 + 0.11375;
    double subx_90 = subx_42*subx_89;
    double subx_91 = states(6,0);
    double subx_92 = 10000.0*subx_91;
    double subx_93 = 920.0*subx_86 + subx_92;
    double subx_94 = subx_87 + subx_91;
    double subx_95 = pow(M_PI, 2);
    double subx_96 = 4302.0*subx_95;
    double subx_97 = 0.11375*subx_72;
    double subx_98 = subx_72*subx_88;
    double subx_99 = subx_97 + subx_98;
    double subx_100 = 0.11375*subx_76;
    double subx_101 = -subx_76*subx_88;
    double subx_102 = -subx_100 + subx_101;
    double subx_103 = 43.02*M_PI;
    double subx_104 = (-subx_103*(subx_102*subx_66 + subx_45*subx_99 + subx_81) - subx_96*(subx_90 + subx_94))*(0.5*erf(10000.0*subx_90 + subx_93) + 0.5);
    double subx_105 = subx_74 + subx_84;
    double subx_106 = subx_35 + subx_36 + subx_37 + subx_44;
    double subx_107 = subx_106*subx_72 + subx_25*subx_76;
    double subx_108 = subx_102*subx_59 + subx_31*subx_99 + subx_80 - subx_87*(subx_105*subx_20 + subx_107);
    double subx_109 = subx_31*subx_76 + subx_59*subx_72;
    double subx_110 = subx_102*subx_106 + subx_25*subx_99 + subx_79 + subx_87*(subx_105*subx_39 + subx_109);
    double subx_111 = sqrt(pow(fabs(subx_108), 2) + pow(fabs(subx_110), 2));
    double subx_112 = 1.0/subx_111;
    double subx_113 = ((subx_111 > 0.0) ? (   subx_108*subx_112): (   0.0));
    double subx_114 = -subx_104;
    double subx_115 = ((subx_114 < 0) ? (   0): (   subx_114));
    double subx_116 = ((subx_111 < 0.0001) ? (   12000.0*subx_111): ((subx_111 >= 0.0001 && subx_111 < 0.0002) ? (   -2000.0*subx_111 + 1.4): (   1.0)));
    double subx_117 = subx_115*subx_116;
    double subx_118 = subx_113*subx_117;
    double subx_119 = ((subx_111 > 0.0) ? (   subx_110*subx_112): (   1.0));
    double subx_120 = subx_117*subx_119;
    double subx_121 = subx_88 - 0.11375;
    double subx_122 = subx_121*subx_42;
    double subx_123 = subx_100 + subx_101;
    double subx_124 = -subx_97 + subx_98;
    double subx_125 = (-subx_103*(subx_123*subx_66 + subx_124*subx_45 + subx_81) - subx_96*(subx_122 + subx_94))*(0.5*erf(10000.0*subx_122 + subx_93) + 0.5);
    double subx_126 = subx_74 + subx_83;
    double subx_127 = subx_123*subx_59 + subx_124*subx_31 + subx_80 - subx_87*(subx_107 + subx_126*subx_20);
    double subx_128 = subx_106*subx_123 + subx_124*subx_25 + subx_79 + subx_87*(subx_109 + subx_126*subx_39);
    double subx_129 = sqrt(pow(fabs(subx_127), 2) + pow(fabs(subx_128), 2));
    double subx_130 = 1.0/subx_129;
    double subx_131 = ((subx_129 > 0.0) ? (   subx_127*subx_130): (   0.0));
    double subx_132 = -subx_125;
    double subx_133 = ((subx_132 < 0) ? (   0): (   subx_132));
    double subx_134 = ((subx_129 < 0.0001) ? (   12000.0*subx_129): ((subx_129 >= 0.0001 && subx_129 < 0.0002) ? (   -2000.0*subx_129 + 1.4): (   1.0)));
    double subx_135 = subx_133*subx_134;
    double subx_136 = subx_131*subx_135;
    double subx_137 = ((subx_129 > 0.0) ? (   subx_128*subx_130): (   1.0));
    double subx_138 = subx_135*subx_137;
    double subx_139 = 0.000804*subx_74;
    double subx_140 = subx_139 + 0.000804*subx_83;
    double subx_141 = subx_1*subx_72;
    double subx_142 = subx_4*subx_76;
    double subx_143 = subx_141 - subx_142;
    double subx_144 = 0.000402*subx_141 - 0.000402*subx_142;
    double subx_145 = 0.000402*subx_126*subx_72 - 0.000402*subx_72*subx_74;
    double subx_146 = -0.000402*subx_126*subx_76 + 0.000402*subx_74*subx_76;
    double subx_147 = subx_1*subx_145 - subx_126*subx_144 + subx_140*subx_143 + subx_146*subx_4;
    double subx_148 = subx_139 + 0.000804*subx_84;
    double subx_149 = subx_3*subx_72;
    double subx_150 = subx_5*subx_76;
    double subx_151 = subx_149 - subx_150;
    double subx_152 = 0.000402*subx_149 - 0.000402*subx_150;
    double subx_153 = subx_72*subx_74;
    double subx_154 = subx_105*subx_72 - subx_153;
    double subx_155 = 0.000402*subx_3;
    double subx_156 = subx_74*subx_76;
    double subx_157 = -subx_105*subx_76 + subx_156;
    double subx_158 = 0.000402*subx_5;
    double subx_159 = -subx_105*subx_152 + subx_148*subx_151 + subx_154*subx_155 + subx_157*subx_158;
    double subx_160 = subx_4*subx_72;
    double subx_161 = subx_1*subx_76;
    double subx_162 = 0.000402*subx_160 + 0.000402*subx_161;
    double subx_163 = subx_160 + subx_161;
    double subx_164 = subx_1*subx_146 + subx_126*subx_162 - subx_140*subx_163 - subx_145*subx_4;
    double subx_165 = subx_5*subx_72;
    double subx_166 = subx_3*subx_76;
    double subx_167 = 0.000402*subx_165 + 0.000402*subx_166;
    double subx_168 = subx_165 + subx_166;
    double subx_169 = subx_105*subx_167 - subx_148*subx_168 - subx_154*subx_158 + subx_155*subx_157;
    double subx_170 = 0.092*subx_113*subx_115*subx_116*subx_86;
    double subx_171 = 0.092*subx_131*subx_133*subx_134*subx_86;
    double subx_172 = 0.092*subx_115*subx_116*subx_119*subx_86;
    double subx_173 = 0.092*subx_133*subx_134*subx_137*subx_86;
    double subx_174 = subx_10*subx_76;
    double subx_175 = subx_74*subx_8;
    double subx_176 = subx_174 - subx_175;
    double subx_177 = subx_10*subx_74;
    double subx_178 = subx_76*subx_8;
    double subx_179 = 0.00417793333333333*subx_177 + 0.00417793333333333*subx_178;
    double subx_180 = 0.00937333333333334*subx_177 + 0.00937333333333334*subx_178;
    double subx_181 = 0.108888707562375*subx_177 + 0.108888707562375*subx_178;
    double subx_182 = subx_177 + subx_178;
    double subx_183 = 9.891512475e-5*subx_174 - 9.891512475e-5*subx_175;
    double subx_184 = 0.00241793333333333*subx_174 - 0.00241793333333333*subx_175;
    double subx_185 = 2.0e-5*subx_174 - 2.0e-5*subx_175;
    double subx_186 = -0.023645*subx_177 - 0.023645*subx_178;
    double subx_187 = subx_72 + subx_82;
    double subx_188 = -0.915*subx_177 - 0.915*subx_178;
    double subx_189 = -0.265*subx_177 - 0.265*subx_178;
    double subx_190 = -1.705*subx_177 - 1.705*subx_178;
    double subx_191 = 0.03225178*subx_82;
    double subx_192 = 0.5*erf(-18300.0*subx_63 + 18300.0*subx_64 + subx_92) + 0.5;
    double subx_193 = -subx_103*(subx_62*(1.83*subx_72 + 1.83*subx_82) + subx_66*(-1.83*subx_177 - 1.83*subx_178) + subx_81) - subx_96*(-1.83*subx_63 + 1.83*subx_64 + subx_91);
    double subx_194 = subx_176*subx_179 + subx_176*subx_180 + subx_176*subx_181 - 0.35685*subx_176*subx_188 - 0.106*subx_176*subx_189 - 1.364*subx_176*subx_190 + subx_177*subx_191 + subx_178*subx_191 - subx_182*subx_183 - subx_182*subx_184 - subx_182*subx_185 + 1.364*subx_186*subx_187 + 1.83*subx_192*subx_193*subx_62 + 17.9150958675*subx_43 + 17.9150958675*subx_46;
    double subx_195 = subx_151*subx_167;
    double subx_196 = 0.00657333333333333*subx_72 + 0.00657333333333333*subx_82;
    double subx_197 = subx_153 - subx_187*subx_74;
    double subx_198 = subx_10*subx_197;
    double subx_199 = subx_187*subx_76 - subx_72*subx_76;
    double subx_200 = subx_199*subx_8;
    double subx_201 = subx_179*subx_187 - subx_182*subx_196 + 0.00241793333333333*subx_198 - 0.00241793333333333*subx_200;
    double subx_202 = 0.00937333333333334*subx_72 + 0.00937333333333334*subx_82;
    double subx_203 = subx_180*subx_187 - subx_182*subx_202 + 2.0e-5*subx_198 - 2.0e-5*subx_200;
    double subx_204 = 0.108888707562375*subx_72 + 0.108888707562375*subx_82;
    double subx_205 = subx_181*subx_187 - subx_182*subx_204 + 9.891512475e-5*subx_198 - 9.891512475e-5*subx_200;
    double subx_206 = -subx_152*subx_168;
    double subx_207 = subx_10*subx_199;
    double subx_208 = subx_197*subx_8;
    double subx_209 = subx_176*subx_202 - subx_185*subx_187 + 0.00937333333333334*subx_207 + 0.00937333333333334*subx_208;
    double subx_210 = subx_176*subx_196 - subx_184*subx_187 + 0.00417793333333333*subx_207 + 0.00417793333333333*subx_208;
    double subx_211 = subx_176*subx_204 - subx_183*subx_187 + 0.108888707562375*subx_207 + 0.108888707562375*subx_208;
    double subx_212 = 0.8*subx_10*subx_76 - 0.8*subx_74*subx_8;
    double subx_213 = -0.018916*subx_10*subx_74 - 0.018916*subx_76*subx_8;
    double subx_214 = 0.018916*subx_82;
    double subx_215 = -subx_177*subx_214 - subx_178*subx_214 - subx_187*subx_213 + subx_190*subx_212;
    double subx_216 = 0.023645*subx_8;
    double subx_217 = 0.023645*subx_174 - 0.023645*subx_175 + 1.705*subx_72 + 1.705*subx_82;
    double subx_218 = -subx_174*subx_214 + subx_175*subx_214 - 0.8*subx_182*subx_190 + 0.8*subx_187*subx_217;
    double subx_219 = 0.023645*subx_10;
    double subx_220 = -2.32562*subx_10*subx_76*subx_82 - 1.364*subx_176*subx_217 + 1.364*subx_182*subx_186 + 2.32562*subx_74*subx_8*subx_82;
    double subx_221 = 0.35685*subx_72 + 0.35685*subx_82;
    double subx_222 = subx_74*subx_8*subx_82;
    double subx_223 = subx_10*subx_76*subx_82;
    double subx_224 = -subx_176*subx_221 + 0.35685*subx_222 - 0.35685*subx_223;
    double subx_225 = 0.915*subx_10;
    double subx_226 = -0.02809*subx_10*subx_76*subx_82 - 0.106*subx_176*(0.265*subx_72 + 0.265*subx_82) + 0.02809*subx_74*subx_8*subx_82;
    double subx_227 = 7.84524*subx_43 + 7.84524*subx_46;
    double subx_228 = 7.84524*subx_63 - 7.84524*subx_64;
    double subx_229 = 26.7522684*subx_12*subx_16 - 26.7522684*subx_13*subx_17;
    double subx_230 = -7.649109*subx_13*subx_17 + 7.649109*subx_23;
    double subx_231 = 2.0789886*subx_12*subx_16 - 2.0789886*subx_13*subx_17;
    double subx_232 = 1.83*subx_192*subx_193*subx_66;
    double subx_233 = subx_170*subx_20;
    double subx_234 = -subx_172*subx_39;
    double subx_235 = subx_143*subx_162 - subx_144*subx_163 + subx_171*subx_20 - subx_173*subx_39;
    double subx_236 = -subx_88;
    double subx_237 = 0.915*subx_8;
    double subx_238 = 0.092*subx_25*subx_86;
    double subx_239 = 0.02275*pow(subx_72, 2);
    double subx_240 = 0.02275*pow(subx_76, 2);
    double subx_241 = subx_239 + subx_240;
    double subx_242 = -subx_239 - subx_240;
    double subx_243 = -0.39*subx_182*subx_188 + subx_187*subx_221;
    double subx_244 = 0.106*subx_72 + 0.106*subx_82;
    double subx_245 = -0.4*subx_182*subx_189 + subx_187*subx_244;
    double subx_246 = -subx_176*subx_244 + 0.106*subx_222 - 0.106*subx_223;
    double subx_247 = subx_182*subx_213 - subx_212*subx_217 + 1.364*subx_222 - 1.364*subx_223;
    double subx_248 = 0.39*subx_176*subx_188;
    double subx_249 = 0.4*subx_176*subx_189;

    Eigen::Matrix<double,19,19> ret;
    ret(0,0) = 1;
    ret(0,1) = 0;
    ret(0,2) = 0;
    ret(0,3) = 0;
    ret(0,4) = 0;
    ret(0,5) = 0;
    ret(0,6) = 0;
    ret(0,7) = 0;
    ret(0,8) = 0;
    ret(0,9) = 0;
    ret(0,10) = 0;
    ret(0,11) = 0;
    ret(0,12) = 0;
    ret(0,13) = 0;
    ret(0,14) = 0;
    ret(0,15) = 0;
    ret(0,16) = 0;
    ret(0,17) = 0;
    ret(0,18) = 0;
    ret(1,0) = 0;
    ret(1,1) = 1;
    ret(1,2) = 0;
    ret(1,3) = 0;
    ret(1,4) = 0;
    ret(1,5) = 0;
    ret(1,6) = 0;
    ret(1,7) = 0;
    ret(1,8) = 0;
    ret(1,9) = 0;
    ret(1,10) = 0;
    ret(1,11) = 0;
    ret(1,12) = 0;
    ret(1,13) = 0;
    ret(1,14) = 0;
    ret(1,15) = 0;
    ret(1,16) = 0;
    ret(1,17) = 0;
    ret(1,18) = 0;
    ret(2,0) = 0;
    ret(2,1) = 0;
    ret(2,2) = 1;
    ret(2,3) = 0;
    ret(2,4) = 0;
    ret(2,5) = 0;
    ret(2,6) = 0;
    ret(2,7) = 0;
    ret(2,8) = 0;
    ret(2,9) = 0;
    ret(2,10) = 0;
    ret(2,11) = 0;
    ret(2,12) = 0;
    ret(2,13) = 0;
    ret(2,14) = 0;
    ret(2,15) = 0;
    ret(2,16) = 0;
    ret(2,17) = 0;
    ret(2,18) = 0;
    ret(3,0) = 0;
    ret(3,1) = 0;
    ret(3,2) = 0;
    ret(3,3) = 1;
    ret(3,4) = 0;
    ret(3,5) = 0;
    ret(3,6) = 0;
    ret(3,7) = 0;
    ret(3,8) = 0;
    ret(3,9) = 0;
    ret(3,10) = 0;
    ret(3,11) = 0;
    ret(3,12) = 0;
    ret(3,13) = 0;
    ret(3,14) = 0;
    ret(3,15) = 0;
    ret(3,16) = 0;
    ret(3,17) = 0;
    ret(3,18) = 0;
    ret(4,0) = 0;
    ret(4,1) = 0;
    ret(4,2) = 0;
    ret(4,3) = 0;
    ret(4,4) = 1;
    ret(4,5) = 0;
    ret(4,6) = 0;
    ret(4,7) = 0;
    ret(4,8) = 0;
    ret(4,9) = 0;
    ret(4,10) = 0;
    ret(4,11) = 0;
    ret(4,12) = 0;
    ret(4,13) = 0;
    ret(4,14) = 0;
    ret(4,15) = 0;
    ret(4,16) = 0;
    ret(4,17) = 0;
    ret(4,18) = 0;
    ret(5,0) = 0;
    ret(5,1) = 0;
    ret(5,2) = 0;
    ret(5,3) = 0;
    ret(5,4) = 0;
    ret(5,5) = 1;
    ret(5,6) = 0;
    ret(5,7) = 0;
    ret(5,8) = 0;
    ret(5,9) = 0;
    ret(5,10) = 0;
    ret(5,11) = 0;
    ret(5,12) = 0;
    ret(5,13) = 0;
    ret(5,14) = 0;
    ret(5,15) = 0;
    ret(5,16) = 0;
    ret(5,17) = 0;
    ret(5,18) = 0;
    ret(6,0) = 0;
    ret(6,1) = 0;
    ret(6,2) = 0;
    ret(6,3) = 0;
    ret(6,4) = 0;
    ret(6,5) = 0;
    ret(6,6) = 1;
    ret(6,7) = 0;
    ret(6,8) = 0;
    ret(6,9) = 0;
    ret(6,10) = 0;
    ret(6,11) = 0;
    ret(6,12) = 0;
    ret(6,13) = 0;
    ret(6,14) = 0;
    ret(6,15) = 0;
    ret(6,16) = 0;
    ret(6,17) = 0;
    ret(6,18) = 0;
    ret(7,0) = 0;
    ret(7,1) = 0;
    ret(7,2) = 0;
    ret(7,3) = 0;
    ret(7,4) = 0;
    ret(7,5) = 0;
    ret(7,6) = 0;
    ret(7,7) = 1;
    ret(7,8) = 0;
    ret(7,9) = 0;
    ret(7,10) = 0;
    ret(7,11) = 0;
    ret(7,12) = 0;
    ret(7,13) = 0;
    ret(7,14) = 0;
    ret(7,15) = 0;
    ret(7,16) = 0;
    ret(7,17) = 0;
    ret(7,18) = 0;
    ret(8,0) = 0;
    ret(8,1) = 0;
    ret(8,2) = 0;
    ret(8,3) = 0;
    ret(8,4) = 0;
    ret(8,5) = 0;
    ret(8,6) = 0;
    ret(8,7) = 0;
    ret(8,8) = 1;
    ret(8,9) = 0;
    ret(8,10) = 0;
    ret(8,11) = 0;
    ret(8,12) = 0;
    ret(8,13) = 0;
    ret(8,14) = 0;
    ret(8,15) = 0;
    ret(8,16) = 0;
    ret(8,17) = 0;
    ret(8,18) = 0;
    ret(9,0) = 0;
    ret(9,1) = 0;
    ret(9,2) = 0;
    ret(9,3) = 0;
    ret(9,4) = 0;
    ret(9,5) = 0;
    ret(9,6) = 0;
    ret(9,7) = 0;
    ret(9,8) = 0;
    ret(9,9) = 1;
    ret(9,10) = 0;
    ret(9,11) = 0;
    ret(9,12) = 0;
    ret(9,13) = 0;
    ret(9,14) = 0;
    ret(9,15) = 0;
    ret(9,16) = 0;
    ret(9,17) = 0;
    ret(9,18) = 0;
    ret(10,0) = 0;
    ret(10,1) = 0;
    ret(10,2) = 0;
    ret(10,3) = 0;
    ret(10,4) = 0;
    ret(10,5) = 0;
    ret(10,6) = 0;
    ret(10,7) = 0;
    ret(10,8) = 0;
    ret(10,9) = 0;
    ret(10,10) = subx_6 + 2.81223888256237;
    ret(10,11) = subx_9;
    ret(10,12) = subx_11;
    ret(10,13) = subx_27;
    ret(10,14) = subx_41;
    ret(10,15) = subx_47;
    ret(10,16) = 2.80506312422904;
    ret(10,17) = 0;
    ret(10,18) = 0;
    ret(11,0) = 0;
    ret(11,1) = 0;
    ret(11,2) = 0;
    ret(11,3) = 0;
    ret(11,4) = 0;
    ret(11,5) = 0;
    ret(11,6) = 0;
    ret(11,7) = 0;
    ret(11,8) = 0;
    ret(11,9) = 0;
    ret(11,10) = subx_9;
    ret(11,11) = 2.80311499304904*subx_48 + 0.00298411727808333*subx_49 + 0.002888;
    ret(11,12) = subx_50;
    ret(11,13) = subx_56;
    ret(11,14) = subx_61;
    ret(11,15) = subx_67;
    ret(11,16) = subx_9;
    ret(11,17) = 0.000804000000000000;
    ret(11,18) = 0.000804000000000000;
    ret(12,0) = 0;
    ret(12,1) = 0;
    ret(12,2) = 0;
    ret(12,3) = 0;
    ret(12,4) = 0;
    ret(12,5) = 0;
    ret(12,6) = 0;
    ret(12,7) = 0;
    ret(12,8) = 0;
    ret(12,9) = 0;
    ret(12,10) = subx_11;
    ret(12,11) = subx_50;
    ret(12,12) = 0.00298411727808333*subx_48 + 2.80311499304904*subx_49 + subx_6 + 0.00717575833333333;
    ret(12,13) = subx_68;
    ret(12,14) = subx_70;
    ret(12,15) = subx_71;
    ret(12,16) = subx_11;
    ret(12,17) = 0;
    ret(12,18) = 0;
    ret(13,0) = 0;
    ret(13,1) = 0;
    ret(13,2) = 0;
    ret(13,3) = 0;
    ret(13,4) = 0;
    ret(13,5) = 0;
    ret(13,6) = 0;
    ret(13,7) = 0;
    ret(13,8) = 0;
    ret(13,9) = 0;
    ret(13,10) = subx_27;
    ret(13,11) = subx_56;
    ret(13,12) = subx_68;
    ret(13,13) = 2.39000000000000;
    ret(13,14) = 0;
    ret(13,15) = 0;
    ret(13,16) = subx_27;
    ret(13,17) = 0;
    ret(13,18) = 0;
    ret(14,0) = 0;
    ret(14,1) = 0;
    ret(14,2) = 0;
    ret(14,3) = 0;
    ret(14,4) = 0;
    ret(14,5) = 0;
    ret(14,6) = 0;
    ret(14,7) = 0;
    ret(14,8) = 0;
    ret(14,9) = 0;
    ret(14,10) = subx_41;
    ret(14,11) = subx_61;
    ret(14,12) = subx_70;
    ret(14,13) = 0;
    ret(14,14) = 2.39000000000000;
    ret(14,15) = 0;
    ret(14,16) = subx_41;
    ret(14,17) = 0;
    ret(14,18) = 0;
    ret(15,0) = 0;
    ret(15,1) = 0;
    ret(15,2) = 0;
    ret(15,3) = 0;
    ret(15,4) = 0;
    ret(15,5) = 0;
    ret(15,6) = 0;
    ret(15,7) = 0;
    ret(15,8) = 0;
    ret(15,9) = 0;
    ret(15,10) = subx_47;
    ret(15,11) = subx_67;
    ret(15,12) = subx_71;
    ret(15,13) = 0;
    ret(15,14) = 0;
    ret(15,15) = 2.39000000000000;
    ret(15,16) = subx_47;
    ret(15,17) = 0;
    ret(15,18) = 0;
    ret(16,0) = 0;
    ret(16,1) = 0;
    ret(16,2) = 0;
    ret(16,3) = 0;
    ret(16,4) = 0;
    ret(16,5) = 0;
    ret(16,6) = 0;
    ret(16,7) = 0;
    ret(16,8) = 0;
    ret(16,9) = 0;
    ret(16,10) = 2.80506312422904;
    ret(16,11) = subx_9;
    ret(16,12) = subx_11;
    ret(16,13) = subx_27;
    ret(16,14) = subx_41;
    ret(16,15) = subx_47;
    ret(16,16) = 2.80506312422904;
    ret(16,17) = 0;
    ret(16,18) = 0;
    ret(17,0) = 0;
    ret(17,1) = 0;
    ret(17,2) = 0;
    ret(17,3) = 0;
    ret(17,4) = 0;
    ret(17,5) = 0;
    ret(17,6) = 0;
    ret(17,7) = 0;
    ret(17,8) = 0;
    ret(17,9) = 0;
    ret(17,10) = 0;
    ret(17,11) = 0.000804000000000000;
    ret(17,12) = 0;
    ret(17,13) = 0;
    ret(17,14) = 0;
    ret(17,15) = 0;
    ret(17,16) = 0;
    ret(17,17) = 0.000804000000000000;
    ret(17,18) = 0;
    ret(18,0) = 0;
    ret(18,1) = 0;
    ret(18,2) = 0;
    ret(18,3) = 0;
    ret(18,4) = 0;
    ret(18,5) = 0;
    ret(18,6) = 0;
    ret(18,7) = 0;
    ret(18,8) = 0;
    ret(18,9) = 0;
    ret(18,10) = 0;
    ret(18,11) = 0.000804000000000000;
    ret(18,12) = 0;
    ret(18,13) = 0;
    ret(18,14) = 0;
    ret(18,15) = 0;
    ret(18,16) = 0;
    ret(18,17) = 0;
    ret(18,18) = 0.000804000000000000;
    return ret;
}

Eigen::Matrix<double,19,1> get_fo(Eigen::Matrix<double,19,1> states, Eigen::Matrix<double,2,1> inputs) {
    double subx_0 = states(8,0);
    double subx_1 = cos(subx_0);
    double subx_2 = states(9,0);
    double subx_3 = cos(subx_2);
    double subx_4 = sin(subx_0);
    double subx_5 = sin(subx_2);
    double subx_6 = 0.000402*pow(subx_1, 2) + 0.000402*pow(subx_3, 2) + 0.000402*pow(subx_4, 2) + 0.000402*pow(subx_5, 2);
    double subx_7 = states(7,0);
    double subx_8 = sin(subx_7);
    double subx_9 = -0.03225178*subx_8;
    double subx_10 = cos(subx_7);
    double subx_11 = 0.03225178*subx_10;
    double subx_12 = states(3,0);
    double subx_13 = states(0,0);
    double subx_14 = 2*subx_13;
    double subx_15 = subx_12*subx_14;
    double subx_16 = states(1,0);
    double subx_17 = states(2,0);
    double subx_18 = 2*subx_17;
    double subx_19 = subx_16*subx_18;
    double subx_20 = -subx_15 + subx_19;
    double subx_21 = subx_10*subx_20;
    double subx_22 = subx_13*subx_18;
    double subx_23 = subx_12*subx_16;
    double subx_24 = 2*subx_23;
    double subx_25 = subx_22 + subx_24;
    double subx_26 = subx_25*subx_8;
    double subx_27 = 1.82685*subx_21 + 1.82685*subx_26;
    double subx_28 = subx_14*subx_16;
    double subx_29 = -subx_28;
    double subx_30 = subx_12*subx_18;
    double subx_31 = subx_29 + subx_30;
    double subx_32 = subx_31*subx_8;
    double subx_33 = pow(subx_17, 2);
    double subx_34 = pow(subx_12, 2);
    double subx_35 = -subx_34;
    double subx_36 = pow(subx_13, 2);
    double subx_37 = pow(subx_16, 2);
    double subx_38 = subx_36 - subx_37;
    double subx_39 = subx_33 + subx_35 + subx_38;
    double subx_40 = subx_10*subx_39;
    double subx_41 = 1.82685*subx_32 + 1.82685*subx_40;
    double subx_42 = subx_28 + subx_30;
    double subx_43 = subx_10*subx_42;
    double subx_44 = -subx_33;
    double subx_45 = subx_34 + subx_38 + subx_44;
    double subx_46 = subx_45*subx_8;
    double subx_47 = 1.82685*subx_43 + 1.82685*subx_46;
    double subx_48 = pow(subx_10, 2);
    double subx_49 = pow(subx_8, 2);
    double subx_50 = 2.80013087577096*subx_10*subx_8;
    double subx_51 = subx_21 + subx_26;
    double subx_52 = 0.018916*subx_8;
    double subx_53 = subx_10*subx_25 - subx_20*subx_8;
    double subx_54 = 0.018916*subx_10;
    double subx_55 = -1.82685*subx_33 - 1.82685*subx_34 + 1.82685*subx_36 + 1.82685*subx_37;
    double subx_56 = -subx_10*subx_55 - subx_51*subx_52 - subx_53*subx_54;
    double subx_57 = subx_10*subx_31 - subx_39*subx_8;
    double subx_58 = subx_32 + subx_40;
    double subx_59 = subx_15 + subx_19;
    double subx_60 = 1.82685*subx_10;
    double subx_61 = -subx_52*subx_58 - subx_54*subx_57 - subx_59*subx_60;
    double subx_62 = subx_43 + subx_46;
    double subx_63 = subx_10*subx_45;
    double subx_64 = subx_42*subx_8;
    double subx_65 = subx_63 - subx_64;
    double subx_66 = -subx_22 + subx_24;
    double subx_67 = -subx_52*subx_62 - subx_54*subx_65 - subx_60*subx_66;
    double subx_68 = subx_51*subx_54 - subx_52*subx_53 - subx_55*subx_8;
    double subx_69 = 1.82685*subx_8;
    double subx_70 = -subx_52*subx_57 + subx_54*subx_58 - subx_59*subx_69;
    double subx_71 = -subx_52*subx_65 + subx_54*subx_62 - subx_66*subx_69;
    double subx_72 = states(10,0);
    double subx_73 = 0.5*subx_72;
    double subx_74 = states(11,0);
    double subx_75 = 0.5*subx_17;
    double subx_76 = states(12,0);
    double subx_77 = 0.5*subx_76;
    double subx_78 = 0.5*subx_74;
    double subx_79 = states(13,0);
    double subx_80 = states(14,0);
    double subx_81 = states(15,0);
    double subx_82 = states(16,0);
    double subx_83 = states(17,0);
    double subx_84 = states(18,0);
    double subx_85 = subx_29 - subx_30;
    double subx_86 = pow(2*subx_42*subx_85 + pow(subx_85, 2) + 1, -1.0L/2.0L);
    double subx_87 = 0.092*subx_86;
    double subx_88 = subx_85*subx_87;
    double subx_89 = subx_88 + 0.11375;
    double subx_90 = subx_42*subx_89;
    double subx_91 = states(6,0);
    double subx_92 = 10000.0*subx_91;
    double subx_93 = 920.0*subx_86 + subx_92;
    double subx_94 = subx_87 + subx_91;
    double subx_95 = pow(M_PI, 2);
    double subx_96 = 4302.0*subx_95;
    double subx_97 = 0.11375*subx_72;
    double subx_98 = subx_72*subx_88;
    double subx_99 = subx_97 + subx_98;
    double subx_100 = 0.11375*subx_76;
    double subx_101 = -subx_76*subx_88;
    double subx_102 = -subx_100 + subx_101;
    double subx_103 = 43.02*M_PI;
    double subx_104 = (-subx_103*(subx_102*subx_66 + subx_45*subx_99 + subx_81) - subx_96*(subx_90 + subx_94))*(0.5*erf(10000.0*subx_90 + subx_93) + 0.5);
    double subx_105 = subx_74 + subx_84;
    double subx_106 = subx_35 + subx_36 + subx_37 + subx_44;
    double subx_107 = subx_106*subx_72 + subx_25*subx_76;
    double subx_108 = subx_102*subx_59 + subx_31*subx_99 + subx_80 - subx_87*(subx_105*subx_20 + subx_107);
    double subx_109 = subx_31*subx_76 + subx_59*subx_72;
    double subx_110 = subx_102*subx_106 + subx_25*subx_99 + subx_79 + subx_87*(subx_105*subx_39 + subx_109);
    double subx_111 = sqrt(pow(fabs(subx_108), 2) + pow(fabs(subx_110), 2));
    double subx_112 = 1.0/subx_111;
    double subx_113 = ((subx_111 > 0.0) ? (   subx_108*subx_112): (   0.0));
    double subx_114 = -subx_104;
    double subx_115 = ((subx_114 < 0) ? (   0): (   subx_114));
    double subx_116 = ((subx_111 < 0.0001) ? (   12000.0*subx_111): ((subx_111 >= 0.0001 && subx_111 < 0.0002) ? (   -2000.0*subx_111 + 1.4): (   1.0)));
    double subx_117 = subx_115*subx_116;
    double subx_118 = subx_113*subx_117;
    double subx_119 = ((subx_111 > 0.0) ? (   subx_110*subx_112): (   1.0));
    double subx_120 = subx_117*subx_119;
    double subx_121 = subx_88 - 0.11375;
    double subx_122 = subx_121*subx_42;
    double subx_123 = subx_100 + subx_101;
    double subx_124 = -subx_97 + subx_98;
    double subx_125 = (-subx_103*(subx_123*subx_66 + subx_124*subx_45 + subx_81) - subx_96*(subx_122 + subx_94))*(0.5*erf(10000.0*subx_122 + subx_93) + 0.5);
    double subx_126 = subx_74 + subx_83;
    double subx_127 = subx_123*subx_59 + subx_124*subx_31 + subx_80 - subx_87*(subx_107 + subx_126*subx_20);
    double subx_128 = subx_106*subx_123 + subx_124*subx_25 + subx_79 + subx_87*(subx_109 + subx_126*subx_39);
    double subx_129 = sqrt(pow(fabs(subx_127), 2) + pow(fabs(subx_128), 2));
    double subx_130 = 1.0/subx_129;
    double subx_131 = ((subx_129 > 0.0) ? (   subx_127*subx_130): (   0.0));
    double subx_132 = -subx_125;
    double subx_133 = ((subx_132 < 0) ? (   0): (   subx_132));
    double subx_134 = ((subx_129 < 0.0001) ? (   12000.0*subx_129): ((subx_129 >= 0.0001 && subx_129 < 0.0002) ? (   -2000.0*subx_129 + 1.4): (   1.0)));
    double subx_135 = subx_133*subx_134;
    double subx_136 = subx_131*subx_135;
    double subx_137 = ((subx_129 > 0.0) ? (   subx_128*subx_130): (   1.0));
    double subx_138 = subx_135*subx_137;
    double subx_139 = 0.000804*subx_74;
    double subx_140 = subx_139 + 0.000804*subx_83;
    double subx_141 = subx_1*subx_72;
    double subx_142 = subx_4*subx_76;
    double subx_143 = subx_141 - subx_142;
    double subx_144 = 0.000402*subx_141 - 0.000402*subx_142;
    double subx_145 = 0.000402*subx_126*subx_72 - 0.000402*subx_72*subx_74;
    double subx_146 = -0.000402*subx_126*subx_76 + 0.000402*subx_74*subx_76;
    double subx_147 = subx_1*subx_145 - subx_126*subx_144 + subx_140*subx_143 + subx_146*subx_4;
    double subx_148 = subx_139 + 0.000804*subx_84;
    double subx_149 = subx_3*subx_72;
    double subx_150 = subx_5*subx_76;
    double subx_151 = subx_149 - subx_150;
    double subx_152 = 0.000402*subx_149 - 0.000402*subx_150;
    double subx_153 = subx_72*subx_74;
    double subx_154 = subx_105*subx_72 - subx_153;
    double subx_155 = 0.000402*subx_3;
    double subx_156 = subx_74*subx_76;
    double subx_157 = -subx_105*subx_76 + subx_156;
    double subx_158 = 0.000402*subx_5;
    double subx_159 = -subx_105*subx_152 + subx_148*subx_151 + subx_154*subx_155 + subx_157*subx_158;
    double subx_160 = subx_4*subx_72;
    double subx_161 = subx_1*subx_76;
    double subx_162 = 0.000402*subx_160 + 0.000402*subx_161;
    double subx_163 = subx_160 + subx_161;
    double subx_164 = subx_1*subx_146 + subx_126*subx_162 - subx_140*subx_163 - subx_145*subx_4;
    double subx_165 = subx_5*subx_72;
    double subx_166 = subx_3*subx_76;
    double subx_167 = 0.000402*subx_165 + 0.000402*subx_166;
    double subx_168 = subx_165 + subx_166;
    double subx_169 = subx_105*subx_167 - subx_148*subx_168 - subx_154*subx_158 + subx_155*subx_157;
    double subx_170 = 0.092*subx_113*subx_115*subx_116*subx_86;
    double subx_171 = 0.092*subx_131*subx_133*subx_134*subx_86;
    double subx_172 = 0.092*subx_115*subx_116*subx_119*subx_86;
    double subx_173 = 0.092*subx_133*subx_134*subx_137*subx_86;
    double subx_174 = subx_10*subx_76;
    double subx_175 = subx_74*subx_8;
    double subx_176 = subx_174 - subx_175;
    double subx_177 = subx_10*subx_74;
    double subx_178 = subx_76*subx_8;
    double subx_179 = 0.00417793333333333*subx_177 + 0.00417793333333333*subx_178;
    double subx_180 = 0.00937333333333334*subx_177 + 0.00937333333333334*subx_178;
    double subx_181 = 0.108888707562375*subx_177 + 0.108888707562375*subx_178;
    double subx_182 = subx_177 + subx_178;
    double subx_183 = 9.891512475e-5*subx_174 - 9.891512475e-5*subx_175;
    double subx_184 = 0.00241793333333333*subx_174 - 0.00241793333333333*subx_175;
    double subx_185 = 2.0e-5*subx_174 - 2.0e-5*subx_175;
    double subx_186 = -0.023645*subx_177 - 0.023645*subx_178;
    double subx_187 = subx_72 + subx_82;
    double subx_188 = -0.915*subx_177 - 0.915*subx_178;
    double subx_189 = -0.265*subx_177 - 0.265*subx_178;
    double subx_190 = -1.705*subx_177 - 1.705*subx_178;
    double subx_191 = 0.03225178*subx_82;
    double subx_192 = 0.5*erf(-18300.0*subx_63 + 18300.0*subx_64 + subx_92) + 0.5;
    double subx_193 = -subx_103*(subx_62*(1.83*subx_72 + 1.83*subx_82) + subx_66*(-1.83*subx_177 - 1.83*subx_178) + subx_81) - subx_96*(-1.83*subx_63 + 1.83*subx_64 + subx_91);
    double subx_194 = subx_176*subx_179 + subx_176*subx_180 + subx_176*subx_181 - 0.35685*subx_176*subx_188 - 0.106*subx_176*subx_189 - 1.364*subx_176*subx_190 + subx_177*subx_191 + subx_178*subx_191 - subx_182*subx_183 - subx_182*subx_184 - subx_182*subx_185 + 1.364*subx_186*subx_187 + 1.83*subx_192*subx_193*subx_62 + 17.9150958675*subx_43 + 17.9150958675*subx_46;
    double subx_195 = subx_151*subx_167;
    double subx_196 = 0.00657333333333333*subx_72 + 0.00657333333333333*subx_82;
    double subx_197 = subx_153 - subx_187*subx_74;
    double subx_198 = subx_10*subx_197;
    double subx_199 = subx_187*subx_76 - subx_72*subx_76;
    double subx_200 = subx_199*subx_8;
    double subx_201 = subx_179*subx_187 - subx_182*subx_196 + 0.00241793333333333*subx_198 - 0.00241793333333333*subx_200;
    double subx_202 = 0.00937333333333334*subx_72 + 0.00937333333333334*subx_82;
    double subx_203 = subx_180*subx_187 - subx_182*subx_202 + 2.0e-5*subx_198 - 2.0e-5*subx_200;
    double subx_204 = 0.108888707562375*subx_72 + 0.108888707562375*subx_82;
    double subx_205 = subx_181*subx_187 - subx_182*subx_204 + 9.891512475e-5*subx_198 - 9.891512475e-5*subx_200;
    double subx_206 = -subx_152*subx_168;
    double subx_207 = subx_10*subx_199;
    double subx_208 = subx_197*subx_8;
    double subx_209 = subx_176*subx_202 - subx_185*subx_187 + 0.00937333333333334*subx_207 + 0.00937333333333334*subx_208;
    double subx_210 = subx_176*subx_196 - subx_184*subx_187 + 0.00417793333333333*subx_207 + 0.00417793333333333*subx_208;
    double subx_211 = subx_176*subx_204 - subx_183*subx_187 + 0.108888707562375*subx_207 + 0.108888707562375*subx_208;
    double subx_212 = 0.8*subx_10*subx_76 - 0.8*subx_74*subx_8;
    double subx_213 = -0.018916*subx_10*subx_74 - 0.018916*subx_76*subx_8;
    double subx_214 = 0.018916*subx_82;
    double subx_215 = -subx_177*subx_214 - subx_178*subx_214 - subx_187*subx_213 + subx_190*subx_212;
    double subx_216 = 0.023645*subx_8;
    double subx_217 = 0.023645*subx_174 - 0.023645*subx_175 + 1.705*subx_72 + 1.705*subx_82;
    double subx_218 = -subx_174*subx_214 + subx_175*subx_214 - 0.8*subx_182*subx_190 + 0.8*subx_187*subx_217;
    double subx_219 = 0.023645*subx_10;
    double subx_220 = -2.32562*subx_10*subx_76*subx_82 - 1.364*subx_176*subx_217 + 1.364*subx_182*subx_186 + 2.32562*subx_74*subx_8*subx_82;
    double subx_221 = 0.35685*subx_72 + 0.35685*subx_82;
    double subx_222 = subx_74*subx_8*subx_82;
    double subx_223 = subx_10*subx_76*subx_82;
    double subx_224 = -subx_176*subx_221 + 0.35685*subx_222 - 0.35685*subx_223;
    double subx_225 = 0.915*subx_10;
    double subx_226 = -0.02809*subx_10*subx_76*subx_82 - 0.106*subx_176*(0.265*subx_72 + 0.265*subx_82) + 0.02809*subx_74*subx_8*subx_82;
    double subx_227 = 7.84524*subx_43 + 7.84524*subx_46;
    double subx_228 = 7.84524*subx_63 - 7.84524*subx_64;
    double subx_229 = 26.7522684*subx_12*subx_16 - 26.7522684*subx_13*subx_17;
    double subx_230 = -7.649109*subx_13*subx_17 + 7.649109*subx_23;
    double subx_231 = 2.0789886*subx_12*subx_16 - 2.0789886*subx_13*subx_17;
    double subx_232 = 1.83*subx_192*subx_193*subx_66;
    double subx_233 = subx_170*subx_20;
    double subx_234 = -subx_172*subx_39;
    double subx_235 = subx_143*subx_162 - subx_144*subx_163 + subx_171*subx_20 - subx_173*subx_39;
    double subx_236 = -subx_88;
    double subx_237 = 0.915*subx_8;
    double subx_238 = 0.092*subx_25*subx_86;
    double subx_239 = 0.02275*pow(subx_72, 2);
    double subx_240 = 0.02275*pow(subx_76, 2);
    double subx_241 = subx_239 + subx_240;
    double subx_242 = -subx_239 - subx_240;
    double subx_243 = -0.39*subx_182*subx_188 + subx_187*subx_221;
    double subx_244 = 0.106*subx_72 + 0.106*subx_82;
    double subx_245 = -0.4*subx_182*subx_189 + subx_187*subx_244;
    double subx_246 = -subx_176*subx_244 + 0.106*subx_222 - 0.106*subx_223;
    double subx_247 = subx_182*subx_213 - subx_212*subx_217 + 1.364*subx_222 - 1.364*subx_223;
    double subx_248 = 0.39*subx_176*subx_188;
    double subx_249 = 0.4*subx_176*subx_189;

    Eigen::Matrix<double,19,1> ret;
    ret(0,0) = -subx_12*subx_77 - subx_16*subx_73 - subx_74*subx_75;
    ret(1,0) = -subx_12*subx_78 + subx_13*subx_73 + subx_75*subx_76;
    ret(2,0) = subx_12*subx_73 + subx_13*subx_78 - subx_16*subx_77;
    ret(3,0) = subx_13*subx_77 + subx_16*subx_78 - subx_17*subx_73;
    ret(4,0) = subx_79;
    ret(5,0) = subx_80;
    ret(6,0) = subx_81;
    ret(7,0) = subx_82;
    ret(8,0) = subx_83;
    ret(9,0) = subx_84;
    ret(10,0) = -subx_1*subx_164 + subx_106*subx_170 + subx_106*subx_171 + subx_121*(subx_125*subx_45 - subx_136*subx_31 - subx_138*subx_25) - subx_147*subx_4 - 0.00589575833333333*subx_156 - subx_159*subx_5 - subx_169*subx_3 - subx_172*subx_59 - subx_173*subx_59 + subx_194 + subx_89*(subx_104*subx_45 - subx_118*subx_31 - subx_120*subx_25);
    ret(11,0) = -subx_10*subx_209 - subx_10*subx_210 - subx_10*subx_211 + subx_10*subx_220 + subx_10*subx_226 - subx_10*subx_229 - subx_10*subx_231 - subx_10*subx_232 + subx_195 + subx_201*subx_8 + subx_203*subx_8 + subx_205*subx_8 + subx_206 + subx_215*subx_216 - subx_216*subx_227 + subx_218*subx_219 - subx_219*subx_228 + subx_224*subx_225 - subx_225*subx_230 + subx_233 + subx_234 + subx_235;
    ret(12,0) = -subx_1*subx_147 - subx_10*subx_201 - subx_10*subx_203 - subx_10*subx_205 + subx_118*subx_238 + subx_136*subx_238 + 0.00589575833333333*subx_153 - subx_159*subx_3 + subx_164*subx_4 + subx_169*subx_5 - subx_172*subx_31 - subx_173*subx_31 - subx_209*subx_8 - subx_210*subx_8 - subx_211*subx_8 - subx_215*subx_219 + subx_216*subx_218 - subx_216*subx_228 + subx_219*subx_227 + subx_220*subx_8 + subx_224*subx_237 + subx_226*subx_8 - subx_229*subx_8 - subx_230*subx_237 - subx_231*subx_8 - subx_232*subx_8 + (subx_236 - 0.11375)*(subx_104*subx_66 - subx_106*subx_120 - subx_118*subx_59) + (subx_236 + 0.11375)*(-subx_106*subx_138 + subx_125*subx_66 - subx_136*subx_59);
    ret(13,0) = -subx_106*subx_224 - subx_106*subx_246 - subx_106*subx_247 - subx_120 - subx_138 - subx_20*subx_241 - subx_20*subx_242 - subx_215*subx_51 - subx_218*subx_53 - subx_243*subx_53 - subx_245*subx_53 - subx_248*subx_51 - subx_249*subx_51;
    ret(14,0) = -subx_118 - subx_136 - subx_215*subx_58 - subx_218*subx_57 - subx_224*subx_59 - subx_241*subx_39 - subx_242*subx_39 - subx_243*subx_57 - subx_245*subx_57 - subx_246*subx_59 - subx_247*subx_59 - subx_248*subx_58 - subx_249*subx_58;
    ret(15,0) = subx_104 + subx_125 + subx_192*subx_193 - subx_215*subx_62 - subx_218*subx_65 - subx_224*subx_66 - subx_241*subx_42 - subx_242*subx_42 - subx_243*subx_65 - subx_245*subx_65 - subx_246*subx_66 - subx_247*subx_66 - subx_248*subx_62 - subx_249*subx_62 + 23.4376545;
    ret(16,0) = subx_194 - 25.2455681180614*subx_7*subx_95 - 8.41518937268713*M_PI*subx_82;
    ret(17,0) = inputs(1,0) + subx_235;
    ret(18,0) = inputs(0,0) + subx_195 + subx_206 + subx_233 + subx_234;
    return ret;
}