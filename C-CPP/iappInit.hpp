#include <iostream>
#include <cstring>
template <class actualRealType>
// Iapp values for 100 neurons
void iappInit(std::vector<actualRealType> &iapp, const char *iappOrder) {

    if (iapp.size() != 100)
        return;

    if(!strcmp(iappOrder, "RAND")) {
        // ===== RANDom Order =====
        iapp[0] = to_actualRealType("-2.18619");
        iapp[1] = to_actualRealType("-2.57576");
        iapp[2] = to_actualRealType("-3.41212");
        iapp[3] = to_actualRealType("-3.71471");
        iapp[4] = to_actualRealType("-5.39247");
        iapp[5] = to_actualRealType("2.96518");
        iapp[6] = to_actualRealType("3.7141");
        iapp[7] = to_actualRealType("4.12336");
        iapp[8] = to_actualRealType("2.06839");
        iapp[9] = to_actualRealType("-3.20933");
        iapp[10] = to_actualRealType("-8.78277");
        iapp[11] = to_actualRealType("-8.15378");
        iapp[12] = to_actualRealType("1.66234");
        iapp[13] = to_actualRealType("-5.6827");
        iapp[14] = to_actualRealType("-4.08322");
        iapp[15] = to_actualRealType("-2.69158");
        iapp[16] = to_actualRealType("-2.59407");
        iapp[17] = to_actualRealType("-1.73528");
        iapp[18] = to_actualRealType("2.83197");
        iapp[19] = to_actualRealType("3.40739");
        iapp[20] = to_actualRealType("-5.76464");
        iapp[21] = to_actualRealType("-9.56374");
        iapp[22] = to_actualRealType("-5.77334");
        iapp[23] = to_actualRealType("2.82739");
        iapp[24] = to_actualRealType("-6.84317");
        iapp[25] = to_actualRealType("-6.86834");
        iapp[26] = to_actualRealType("-0.258034");
        iapp[27] = to_actualRealType("-0.471816");
        iapp[28] = to_actualRealType("-3.84289");
        iapp[29] = to_actualRealType("-3.70281");
        iapp[30] = to_actualRealType("-5.67721");
        iapp[31] = to_actualRealType("1.8244");
        iapp[32] = to_actualRealType("-8.66329");
        iapp[33] = to_actualRealType("3.99792");
        iapp[34] = to_actualRealType("-6.82577");
        iapp[35] = to_actualRealType("-5.59847");
        iapp[36] = to_actualRealType("1.23707");
        iapp[37] = to_actualRealType("-3.68862");
        iapp[38] = to_actualRealType("-4.54604");
        iapp[39] = to_actualRealType("4.11283");
        iapp[40] = to_actualRealType("-0.170141");
        iapp[41] = to_actualRealType("0.237281");
        iapp[42] = to_actualRealType("-9.07758");
        iapp[43] = to_actualRealType("-8.10526");
        iapp[44] = to_actualRealType("-8.86425");
        iapp[45] = to_actualRealType("-9.28953");
        iapp[46] = to_actualRealType("-3.72341");
        iapp[47] = to_actualRealType("0.655232");
        iapp[48] = to_actualRealType("-5.05325");
        iapp[49] = to_actualRealType("-4.68108");
        iapp[50] = to_actualRealType("-9.62691");
        iapp[51] = to_actualRealType("-8.53282");
        iapp[52] = to_actualRealType("-7.38884");
        iapp[53] = to_actualRealType("3.81252");
        iapp[54] = to_actualRealType("0.0715659");
        iapp[55] = to_actualRealType("-7.08075");
        iapp[56] = to_actualRealType("-1.63961");
        iapp[57] = to_actualRealType("-3.99533");
        iapp[58] = to_actualRealType("-4.09009");
        iapp[59] = to_actualRealType("-4.03607");
        iapp[60] = to_actualRealType("-5.86123");
        iapp[61] = to_actualRealType("4.15403");
        iapp[62] = to_actualRealType("-0.667745");
        iapp[63] = to_actualRealType("-5.89236");
        iapp[64] = to_actualRealType("-1.58879");
        iapp[65] = to_actualRealType("-1.57781");
        iapp[66] = to_actualRealType("0.690023");
        iapp[67] = to_actualRealType("-8.71548");
        iapp[68] = to_actualRealType("-8.09061");
        iapp[69] = to_actualRealType("-8.4962");
        iapp[70] = to_actualRealType("-0.269478");
        iapp[71] = to_actualRealType("2.26341");
        iapp[72] = to_actualRealType("0.794397");
        iapp[73] = to_actualRealType("1.01733");
        iapp[74] = to_actualRealType("-7.20252");
        iapp[75] = to_actualRealType("2.92398");
        iapp[76] = to_actualRealType("0.197913");
        iapp[77] = to_actualRealType("-4.24802");
        iapp[78] = to_actualRealType("-6.54378");
        iapp[79] = to_actualRealType("-9.74914");
        iapp[80] = to_actualRealType("-1.38554");
        iapp[81] = to_actualRealType("4.37422");
        iapp[82] = to_actualRealType("0.769219");
        iapp[83] = to_actualRealType("4.39207");
        iapp[84] = to_actualRealType("-3.88913");
        iapp[85] = to_actualRealType("2.66762");
        iapp[86] = to_actualRealType("1.93197");
        iapp[87] = to_actualRealType("2.61223");
        iapp[88] = to_actualRealType("-3.17728");
        iapp[89] = to_actualRealType("-2.08136");
        iapp[90] = to_actualRealType("-1.50227");
        iapp[91] = to_actualRealType("-7.49413");
        iapp[92] = to_actualRealType("-6.292");
        iapp[93] = to_actualRealType("1.53645");
        iapp[94] = to_actualRealType("-7.09815");
        iapp[95] = to_actualRealType("-0.705741");
        iapp[96] = to_actualRealType("1.02374");
        iapp[97] = to_actualRealType("-6.51036");
        iapp[98] = to_actualRealType("-0.490585");
        iapp[99] = to_actualRealType("4.55962");
    } else if(!strcmp(iappOrder, "ASC")) {
        //===== ASCeding Order =====
        iapp[0] = to_actualRealType("-9.74914");
        iapp[1] = to_actualRealType("-9.62691");
        iapp[2] = to_actualRealType("-9.56374");
        iapp[3] = to_actualRealType("-9.28953");
        iapp[4] = to_actualRealType("-9.07758");
        iapp[5] = to_actualRealType("-8.86425");
        iapp[6] = to_actualRealType("-8.78277");
        iapp[7] = to_actualRealType("-8.71548");
        iapp[8] = to_actualRealType("-8.66329");
        iapp[9] = to_actualRealType("-8.53282");
        iapp[10] = to_actualRealType("-8.4962");
        iapp[11] = to_actualRealType("-8.15378");
        iapp[12] = to_actualRealType("-8.10526");
        iapp[13] = to_actualRealType("-8.09061");
        iapp[14] = to_actualRealType("-7.49413");
        iapp[15] = to_actualRealType("-7.38884");
        iapp[16] = to_actualRealType("-7.20252");
        iapp[17] = to_actualRealType("-7.09815");
        iapp[18] = to_actualRealType("-7.08075");
        iapp[19] = to_actualRealType("-6.86834");
        iapp[20] = to_actualRealType("-6.84317");
        iapp[21] = to_actualRealType("-6.82577");
        iapp[22] = to_actualRealType("-6.54378");
        iapp[23] = to_actualRealType("-6.51036");
        iapp[24] = to_actualRealType("-6.292");
        iapp[25] = to_actualRealType("-5.89236");
        iapp[26] = to_actualRealType("-5.86123");
        iapp[27] = to_actualRealType("-5.77334");
        iapp[28] = to_actualRealType("-5.76464");
        iapp[29] = to_actualRealType("-5.6827");
        iapp[30] = to_actualRealType("-5.67721");
        iapp[31] = to_actualRealType("-5.59847");
        iapp[32] = to_actualRealType("-5.39247");
        iapp[33] = to_actualRealType("-5.05325");
        iapp[34] = to_actualRealType("-4.68108");
        iapp[35] = to_actualRealType("-4.54604");
        iapp[36] = to_actualRealType("-4.24802");
        iapp[37] = to_actualRealType("-4.09009");
        iapp[38] = to_actualRealType("-4.08322");
        iapp[39] = to_actualRealType("-4.03607");
        iapp[40] = to_actualRealType("-3.99533");
        iapp[41] = to_actualRealType("-3.88913");
        iapp[42] = to_actualRealType("-3.84289");
        iapp[43] = to_actualRealType("-3.72341");
        iapp[44] = to_actualRealType("-3.71471");
        iapp[45] = to_actualRealType("-3.70281");
        iapp[46] = to_actualRealType("-3.68862");
        iapp[47] = to_actualRealType("-3.41212");
        iapp[48] = to_actualRealType("-3.20933");
        iapp[49] = to_actualRealType("-3.17728");
        iapp[50] = to_actualRealType("-2.69158");
        iapp[51] = to_actualRealType("-2.59407");
        iapp[52] = to_actualRealType("-2.57576");
        iapp[53] = to_actualRealType("-2.18619");
        iapp[54] = to_actualRealType("-2.08136");
        iapp[55] = to_actualRealType("-1.73528");
        iapp[56] = to_actualRealType("-1.63961");
        iapp[57] = to_actualRealType("-1.58879");
        iapp[58] = to_actualRealType("-1.57781");
        iapp[59] = to_actualRealType("-1.50227");
        iapp[60] = to_actualRealType("-1.38554");
        iapp[61] = to_actualRealType("-0.705741");
        iapp[62] = to_actualRealType("-0.667745");
        iapp[63] = to_actualRealType("-0.490585");
        iapp[64] = to_actualRealType("-0.471816");
        iapp[65] = to_actualRealType("-0.269478");
        iapp[66] = to_actualRealType("-0.258034");
        iapp[67] = to_actualRealType("-0.170141");
        iapp[68] = to_actualRealType("0.0715659");
        iapp[69] = to_actualRealType("0.197913");
        iapp[70] = to_actualRealType("0.237281");
        iapp[71] = to_actualRealType("0.655232");
        iapp[72] = to_actualRealType("0.690023");
        iapp[73] = to_actualRealType("0.769219");
        iapp[74] = to_actualRealType("0.794397");
        iapp[75] = to_actualRealType("1.01733");
        iapp[76] = to_actualRealType("1.02374");
        iapp[77] = to_actualRealType("1.23707");
        iapp[78] = to_actualRealType("1.53645");
        iapp[79] = to_actualRealType("1.66234");
        iapp[80] = to_actualRealType("1.8244");
        iapp[81] = to_actualRealType("1.93197");
        iapp[82] = to_actualRealType("2.06839");
        iapp[83] = to_actualRealType("2.26341");
        iapp[84] = to_actualRealType("2.61223");
        iapp[85] = to_actualRealType("2.66762");
        iapp[86] = to_actualRealType("2.82739");
        iapp[87] = to_actualRealType("2.83197");
        iapp[88] = to_actualRealType("2.92398");
        iapp[89] = to_actualRealType("2.96518");
        iapp[90] = to_actualRealType("3.40739");
        iapp[91] = to_actualRealType("3.7141");
        iapp[92] = to_actualRealType("3.81252");
        iapp[93] = to_actualRealType("3.99792");
        iapp[94] = to_actualRealType("4.11283");
        iapp[95] = to_actualRealType("4.12336");
        iapp[96] = to_actualRealType("4.15403");
        iapp[97] = to_actualRealType("4.37422");
        iapp[98] = to_actualRealType("4.39207");
        iapp[99] = to_actualRealType("4.55962");
    } else if(!strcmp(iappOrder, "DES")) {
        // ===== DEScending Order =====
        iapp[0] = to_actualRealType("4.55962");
        iapp[1] = to_actualRealType("4.39207");
        iapp[2] = to_actualRealType("4.37422");
        iapp[3] = to_actualRealType("4.15403");
        iapp[4] = to_actualRealType("4.12336");
        iapp[5] = to_actualRealType("4.11283");
        iapp[6] = to_actualRealType("3.99792");
        iapp[7] = to_actualRealType("3.81252");
        iapp[8] = to_actualRealType("3.7141");
        iapp[9] = to_actualRealType("3.40739");
        iapp[10] = to_actualRealType("2.96518");
        iapp[11] = to_actualRealType("2.92398");
        iapp[12] = to_actualRealType("2.83197");
        iapp[13] = to_actualRealType("2.82739");
        iapp[14] = to_actualRealType("2.66762");
        iapp[15] = to_actualRealType("2.61223");
        iapp[16] = to_actualRealType("2.26341");
        iapp[17] = to_actualRealType("2.06839");
        iapp[18] = to_actualRealType("1.93197");
        iapp[19] = to_actualRealType("1.8244");
        iapp[20] = to_actualRealType("1.66234");
        iapp[21] = to_actualRealType("1.53645");
        iapp[22] = to_actualRealType("1.23707");
        iapp[23] = to_actualRealType("1.02374");
        iapp[24] = to_actualRealType("1.01733");
        iapp[25] = to_actualRealType("0.794397");
        iapp[26] = to_actualRealType("0.769219");
        iapp[27] = to_actualRealType("0.690023");
        iapp[28] = to_actualRealType("0.655232");
        iapp[29] = to_actualRealType("0.237281");
        iapp[30] = to_actualRealType("0.197913");
        iapp[31] = to_actualRealType("0.0715659");
        iapp[32] = to_actualRealType("-0.170141");
        iapp[33] = to_actualRealType("-0.258034");
        iapp[34] = to_actualRealType("-0.269478");
        iapp[35] = to_actualRealType("-0.471816");
        iapp[36] = to_actualRealType("-0.490585");
        iapp[37] = to_actualRealType("-0.667745");
        iapp[38] = to_actualRealType("-0.705741");
        iapp[39] = to_actualRealType("-1.38554");
        iapp[40] = to_actualRealType("-1.50227");
        iapp[41] = to_actualRealType("-1.57781");
        iapp[42] = to_actualRealType("-1.58879");
        iapp[43] = to_actualRealType("-1.63961");
        iapp[44] = to_actualRealType("-1.73528");
        iapp[45] = to_actualRealType("-2.08136");
        iapp[46] = to_actualRealType("-2.18619");
        iapp[47] = to_actualRealType("-2.57576");
        iapp[48] = to_actualRealType("-2.59407");
        iapp[49] = to_actualRealType("-2.69158");
        iapp[50] = to_actualRealType("-3.17728");
        iapp[51] = to_actualRealType("-3.20933");
        iapp[52] = to_actualRealType("-3.41212");
        iapp[53] = to_actualRealType("-3.68862");
        iapp[54] = to_actualRealType("-3.70281");
        iapp[55] = to_actualRealType("-3.71471");
        iapp[56] = to_actualRealType("-3.72341");
        iapp[57] = to_actualRealType("-3.84289");
        iapp[58] = to_actualRealType("-3.88913");
        iapp[59] = to_actualRealType("-3.99533");
        iapp[60] = to_actualRealType("-4.03607");
        iapp[61] = to_actualRealType("-4.08322");
        iapp[62] = to_actualRealType("-4.09009");
        iapp[63] = to_actualRealType("-4.24802");
        iapp[64] = to_actualRealType("-4.54604");
        iapp[65] = to_actualRealType("-4.68108");
        iapp[66] = to_actualRealType("-5.05325");
        iapp[67] = to_actualRealType("-5.39247");
        iapp[68] = to_actualRealType("-5.59847");
        iapp[69] = to_actualRealType("-5.67721");
        iapp[70] = to_actualRealType("-5.6827");
        iapp[71] = to_actualRealType("-5.76464");
        iapp[72] = to_actualRealType("-5.77334");
        iapp[73] = to_actualRealType("-5.86123");
        iapp[74] = to_actualRealType("-5.89236");
        iapp[75] = to_actualRealType("-6.292");
        iapp[76] = to_actualRealType("-6.51036");
        iapp[77] = to_actualRealType("-6.54378");
        iapp[78] = to_actualRealType("-6.82577");
        iapp[79] = to_actualRealType("-6.84317");
        iapp[80] = to_actualRealType("-6.86834");
        iapp[81] = to_actualRealType("-7.08075");
        iapp[82] = to_actualRealType("-7.09815");
        iapp[83] = to_actualRealType("-7.20252");
        iapp[84] = to_actualRealType("-7.38884");
        iapp[85] = to_actualRealType("-7.49413");
        iapp[86] = to_actualRealType("-8.09061");
        iapp[87] = to_actualRealType("-8.10526");
        iapp[88] = to_actualRealType("-8.15378");
        iapp[89] = to_actualRealType("-8.4962");
        iapp[90] = to_actualRealType("-8.53282");
        iapp[91] = to_actualRealType("-8.66329");
        iapp[92] = to_actualRealType("-8.71548");
        iapp[93] = to_actualRealType("-8.78277");
        iapp[94] = to_actualRealType("-8.86425");
        iapp[95] = to_actualRealType("-9.07758");
        iapp[96] = to_actualRealType("-9.28953");
        iapp[97] = to_actualRealType("-9.56374");
        iapp[98] = to_actualRealType("-9.62691");
        iapp[99] = to_actualRealType("-9.74914");
    }
}