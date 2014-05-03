#include "../interface/RunHelper.h"

#include "TString.h"





void RunHelper::getBeamPosition( const std::string& runName, float& beamX, float& beamY ) {

  TString runName_tstr(runName.c_str());


  float xTable=-999.;
  float yTable=-999.;
       if( runName_tstr.BeginsWith("BTF_51_" ) ) { xTable = 511.0; yTable = 202.5; }
  else if( runName_tstr.BeginsWith("BTF_52_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_53_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_54_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_55_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_56_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_57_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_58_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_59_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_60_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_61_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_62_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_63_" ) ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_64_" ) ) { xTable = 511.0; yTable = 202.5; }
  else if( runName_tstr.BeginsWith("BTF_65_" ) ) { xTable = 511.0; yTable = 205.0; }
  else if( runName_tstr.BeginsWith("BTF_66_" ) ) { xTable = 511.0; yTable = 197.5; }
  else if( runName_tstr.BeginsWith("BTF_67_" ) ) { xTable = 511.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_68_" ) ) { xTable = 511.0; yTable = 199.0; }
  else if( runName_tstr.BeginsWith("BTF_69_" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_70_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_71_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_72_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_73_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_74_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_75_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_76_" ) ) { xTable = 488.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_77_" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_78_" ) ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_79_" ) ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_80_" ) ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_81_" ) ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_82_" ) ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_83_" ) ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_84_" ) ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_85_" ) ) { xTable = 534.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_86_" ) ) { xTable = 534.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_87_" ) ) { xTable = 534.0; yTable = 222.0; }
  else if( runName_tstr.BeginsWith("BTF_88_" ) ) { xTable = 534.0; yTable = 222.0; }
  else if( runName_tstr.BeginsWith("BTF_89_" ) ) { xTable = 512.0; yTable = 224.0; }
  else if( runName_tstr.BeginsWith("BTF_90_" ) ) { xTable = 490.0; yTable = 224.0; }
  else if( runName_tstr.BeginsWith("BTF_91_" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_92_" ) ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_93_" ) ) { xTable = 514.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_94_" ) ) { xTable = 514.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_95_" ) ) { xTable = 517.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_96_" ) ) { xTable = 517.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_97_" ) ) { xTable = 520.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_98_" ) ) { xTable = 520.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_99_" ) ) { xTable = 520.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_100_") ) { xTable = 520.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_103_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_116_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_117_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_118_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_119_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_120_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_131_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_132_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_133_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_134_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_135_") ) { xTable = 517.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_136_") ) { xTable = 514.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_137_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_138_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_139_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_140_") ) { xTable = 508.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_141_") ) { xTable = 508.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_142_") ) { xTable = 505.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_143_") ) { xTable = 505.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_144_") ) { xTable = 499.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_145_") ) { xTable = 501.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_146_") ) { xTable = 503.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_147_") ) { xTable = 505.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_148_") ) { xTable = 507.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_149_") ) { xTable = 509.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_150_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_151_") ) { xTable = 513.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_152_") ) { xTable = 515.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_153_") ) { xTable = 517.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_154_") ) { xTable = 519.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_155_") ) { xTable = 521.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_156_") ) { xTable = 523.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_157_") ) { xTable = 511.0; yTable = 193.0; }
  else if( runName_tstr.BeginsWith("BTF_158_") ) { xTable = 511.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_159_") ) { xTable = 511.0; yTable = 197.0; }
  else if( runName_tstr.BeginsWith("BTF_160_") ) { xTable = 511.0; yTable = 199.0; }
  else if( runName_tstr.BeginsWith("BTF_161_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_162_") ) { xTable = 511.0; yTable = 203.0; }
  else if( runName_tstr.BeginsWith("BTF_163_") ) { xTable = 511.0; yTable = 205.0; }
  else if( runName_tstr.BeginsWith("BTF_164_") ) { xTable = 511.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_165_") ) { xTable = 511.0; yTable = 209.0; }
  else if( runName_tstr.BeginsWith("BTF_166_") ) { xTable = 505.0; yTable = 209.0; }
  else if( runName_tstr.BeginsWith("BTF_167_") ) { xTable = 502.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_168_") ) { xTable = 502.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_169_") ) { xTable = 502.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_170_") ) { xTable = 502.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_171_") ) { xTable = 502.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_172_") ) { xTable = 505.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_173_") ) { xTable = 505.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_174_") ) { xTable = 505.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_175_") ) { xTable = 508.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_176_") ) { xTable = 508.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_177_") ) { xTable = 508.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_178_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_179_") ) { xTable = 511.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_180_") ) { xTable = 511.0; yTable = 199.0; }
  else if( runName_tstr.BeginsWith("BTF_181_") ) { xTable = 511.0; yTable = 200.0; }
  else if( runName_tstr.BeginsWith("BTF_182_") ) { xTable = 511.0; yTable = 202.0; }
  else if( runName_tstr.BeginsWith("BTF_183_") ) { xTable = 511.0; yTable = 203.0; }
  else if( runName_tstr.BeginsWith("BTF_184_") ) { xTable = 511.0; yTable = 203.0; }
  else if( runName_tstr.BeginsWith("BTF_185_") ) { xTable = 511.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_186_") ) { xTable = 514.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_187_") ) { xTable = 513.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_188_") ) { xTable = 512.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_189_") ) { xTable = 510.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_190_") ) { xTable = 509.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_191_") ) { xTable = 508.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_192_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_193_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_194_") ) { xTable = 511.0; yTable = 186.0; }
  else if( runName_tstr.BeginsWith("BTF_195_") ) { xTable = 511.0; yTable = 189.0; }
  else if( runName_tstr.BeginsWith("BTF_196_") ) { xTable = 511.0; yTable = 192.0; }
  else if( runName_tstr.BeginsWith("BTF_197_") ) { xTable = 511.0; yTable = 195.0; }
  else if( runName_tstr.BeginsWith("BTF_198_") ) { xTable = 511.0; yTable = 198.0; }
  else if( runName_tstr.BeginsWith("BTF_199_") ) { xTable = 511.0; yTable = 204.0; }
  else if( runName_tstr.BeginsWith("BTF_200_") ) { xTable = 511.0; yTable = 207.0; }
  else if( runName_tstr.BeginsWith("BTF_201_") ) { xTable = 511.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_202_") ) { xTable = 511.0; yTable = 210.0; }
  else if( runName_tstr.BeginsWith("BTF_203_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_204_") ) { xTable = 511.0; yTable = 213.0; }
  else if( runName_tstr.BeginsWith("BTF_205_") ) { xTable = 511.0; yTable = 216.0; }
  else if( runName_tstr.BeginsWith("BTF_206_") ) { xTable = 526.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_207_") ) { xTable = 523.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_208_") ) { xTable = 520.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_209_") ) { xTable = 517.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_210_") ) { xTable = 514.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_211_") ) { xTable = 508.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_212_") ) { xTable = 508.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_213_") ) { xTable = 505.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_214_") ) { xTable = 502.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_215_") ) { xTable = 499.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_216_") ) { xTable = 496.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_217_") ) { xTable = 496.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_218_") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_219_") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_220_") ) { xTable = 499.2; yTable = 189.2; }
  else if( runName_tstr.BeginsWith("BTF_221_") ) { xTable = 498.7; yTable = 188.7; }
  else if( runName_tstr.BeginsWith("BTF_222_") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_223_") ) { xTable = 499.7; yTable = 189.7; }
  else if( runName_tstr.BeginsWith("BTF_224_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_225_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_226_") ) { xTable = 496.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_227_") ) { xTable = 511.0; yTable = 201.0; }
  else if( runName_tstr.BeginsWith("BTF_228_") ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_229_") ) { xTable = 532.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_230_") ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_231_") ) { xTable = 510.0; yTable = 178.0; }
  else if( runName_tstr.BeginsWith("BTF_232_") ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_233_") ) { xTable = 488.0; yTable = 180.0; }
  else if( runName_tstr.BeginsWith("BTF_234_") ) { xTable = 488.0; yTable = 202.1; }
  else if( runName_tstr.BeginsWith("BTF_235_") ) { xTable = 488.0; yTable = 202.1; }
  else if( runName_tstr.BeginsWith("BTF_236_") ) { xTable = 534.0; yTable = 200.1; }
  else if( runName_tstr.BeginsWith("BTF_237_") ) { xTable = 534.0; yTable = 200.1; }
  else if( runName_tstr.BeginsWith("BTF_238_") ) { xTable = 534.0; yTable = 222.1; }
  else if( runName_tstr.BeginsWith("BTF_239_") ) { xTable = 534.0; yTable = 222.1; }
  else if( runName_tstr.BeginsWith("BTF_240_") ) { xTable = 512.0; yTable = 224.1; }
  else if( runName_tstr.BeginsWith("BTF_241_") ) { xTable = 512.0; yTable = 224.1; }
  else if( runName_tstr.BeginsWith("BTF_242_") ) { xTable = 512.0; yTable = 224.1; }
  else if( runName_tstr.BeginsWith("BTF_243_") ) { xTable = 490.0; yTable = 224.1; }
  else if( runName_tstr.BeginsWith("BTF_244_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_245_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_246_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_247_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_248_") ) { xTable = 504.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_249_") ) { xTable = 504.0; yTable = 200.6; }
  else if( runName_tstr.BeginsWith("BTF_250_") ) { xTable = 504.0; yTable = 206.6; }
  else if( runName_tstr.BeginsWith("BTF_251_") ) { xTable = 501.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_252_") ) { xTable = 501.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_253_") ) { xTable = 507.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_254_") ) { xTable = 507.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_255_") ) { xTable = 507.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_256_") ) { xTable = 507.0; yTable = 203.6; }
  else if( runName_tstr.BeginsWith("BTF_257_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_258_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_259_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_260_") ) { xTable = 511.0; yTable = 201.1; }
  else if( runName_tstr.BeginsWith("BTF_261_") ) { xTable = 501.2; yTable = 190.2; }
  else if( runName_tstr.BeginsWith("BTF_262_") ) { xTable = 501.2; yTable = 190.2; }
  else if( runName_tstr.BeginsWith("BTF_263_") ) { xTable = 502  ; yTable = 196.5; }
  else if( runName_tstr.BeginsWith("BTF_264_") ) { xTable = 506.5; yTable = 192; }
  else if( runName_tstr.BeginsWith("BTF_265_") ) { xTable = 522.3; yTable = 188.7; }
  else if( runName_tstr.BeginsWith("BTF_266_") ) { xTable = 522.3; yTable = 188.7; }
  else if( runName_tstr.BeginsWith("BTF_267_") ) { xTable = 522.3; yTable = 212.3; }
  else if( runName_tstr.BeginsWith("BTF_268_") ) { xTable = 522.3; yTable = 212.3; }
  else if( runName_tstr.BeginsWith("BTF_269_") ) { xTable = 499.7; yTable = 212.3; }
  else if( runName_tstr.BeginsWith("BTF_270_") ) { xTable = 499.7; yTable = 212.3; }
  else if( runName_tstr.BeginsWith("BTF_271_") ) { xTable = 499.2; yTable = 211.8; }
  else if( runName_tstr.BeginsWith("BTF_272_") ) { xTable = 499.2; yTable = 211.8; }
  else if( runName_tstr.BeginsWith("BTF_273_") ) { xTable = 499.2; yTable = 212.8; }
  else if( runName_tstr.BeginsWith("BTF_274_") ) { xTable = 499.2; yTable = 212.8; }
  else if( runName_tstr.BeginsWith("BTF_275_") ) { xTable = 500.2; yTable = 212.8; }
  else if( runName_tstr.BeginsWith("BTF_276_") ) { xTable = 500.2; yTable = 212.8; }
  else if( runName_tstr.BeginsWith("BTF_277_") ) { xTable = 500.2; yTable = 211.8; }
  else if( runName_tstr.BeginsWith("BTF_278_") ) { xTable = 500.2; yTable = 211.8; }
  else if( runName_tstr.BeginsWith("BTF_279_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_280_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_281_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_282_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_283_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_284_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_285_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_286_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_287_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_288_") ) { xTable = 514.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_289_") ) { xTable = 508.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_290_") ) { xTable = 511.0; yTable = 198; }
  else if( runName_tstr.BeginsWith("BTF_291_") ) { xTable = 511.0; yTable = 204; }
  else if( runName_tstr.BeginsWith("BTF_292_") ) { xTable = 511.0; yTable = 201; }
  else if( runName_tstr.BeginsWith("BTF_293_") ) { xTable = 511.0; yTable = 201; }


  beamX = -xTable + 511.;
  beamY = -yTable + 201.;

}



void RunHelper::getBGOCoordinates( int iChannel, float& x , float& y ) {

  if( iChannel==0 ) {
    x = -532.0+511.; y = -178.0+201.;
  } else if( iChannel==1 ) {
    x = -510.0+511.; y = -178.0+201.;
  } else if( iChannel==2 ) {
    x = -488.0+511.; y = -180.0+201.;
  } else if( iChannel==3 ) {
    x = -534.0+511.; y = -200.1+201.;
  } else if( iChannel==4 ) {
    x = -488.0+511.; y = -202.1+201.;
  } else if( iChannel==5 ) {
    x = -534.0+511.; y = -222.1+201.;
  } else if( iChannel==6 ) {
    x = -512.0+511.; y = -224.1+201.;
  } else if( iChannel==7 ) {
    x = -490.0+511.; y = -224.1+201.;
  }

}
