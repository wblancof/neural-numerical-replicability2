

function iapp = iappInit( iappOrder )

  iapp = zeros(1,100);
  % Random Sorted
  if iappOrder == "RAND" 
    iapp(1) = -2.18619;
    iapp(2) = -2.57576;
    iapp(3) = -3.41212;
    iapp(4) = -3.71471;
    iapp(5) = -5.39247;
    iapp(6) = 2.96518;
    iapp(7) = 3.7141;
    iapp(8) = 4.12336;
    iapp(9) = 2.06839;
    iapp(10) = -3.20933;
    iapp(11) = -8.78277;
    iapp(12) = -8.15378;
    iapp(13) = 1.66234;
    iapp(14) = -5.6827;
    iapp(15) = -4.08322;
    iapp(16) = -2.69158;
    iapp(17) = -2.59407;
    iapp(18) = -1.73528;
    iapp(19) = 2.83197;
    iapp(20) = 3.40739;
    iapp(21) = -5.76464;
    iapp(22) = -9.56374;
    iapp(23) = -5.77334;
    iapp(24) = 2.82739;
    iapp(25) = -6.84317;
    iapp(26) = -6.86834;
    iapp(27) = -0.258034;
    iapp(28) = -0.471816;
    iapp(29) = -3.84289;
    iapp(30) = -3.70281;
    iapp(31) = -5.67721;
    iapp(32) = 1.8244;
    iapp(33) = -8.66329;
    iapp(34) = 3.99792;
    iapp(35) = -6.82577;
    iapp(36) = -5.59847;
    iapp(37) = 1.23707;
    iapp(38) = -3.68862;
    iapp(39) = -4.54604;
    iapp(40) = 4.11283;
    iapp(41) = -0.170141;
    iapp(42) = 0.237281;
    iapp(43) = -9.07758;
    iapp(44) = -8.10526;
    iapp(45) = -8.86425;
    iapp(46) = -9.28953;
    iapp(47) = -3.72341;
    iapp(48) = 0.655232;
    iapp(49) = -5.05325;
    iapp(50) = -4.68108;
    iapp(51) = -9.62691;
    iapp(52) = -8.53282;
    iapp(53) = -7.38884;
    iapp(54) = 3.81252;
    iapp(55) = 0.0715659;
    iapp(56) = -7.08075;
    iapp(57) = -1.63961;
    iapp(58) = -3.99533;
    iapp(59) = -4.09009;
    iapp(60) = -4.03607;
    iapp(61) = -5.86123;
    iapp(62) = 4.15403;
    iapp(63) = -0.667745;
    iapp(64) = -5.89236;
    iapp(65) = -1.58879;
    iapp(66) = -1.57781;
    iapp(67) = 0.690023;
    iapp(68) = -8.71548;
    iapp(69) = -8.09061;
    iapp(70) = -8.4962;
    iapp(71) = -0.269478;
    iapp(72) = 2.26341;
    iapp(73) = 0.794397;
    iapp(74) = 1.01733;
    iapp(75) = -7.20252;
    iapp(76) = 2.92398;
    iapp(77) = 0.197913;
    iapp(78) = -4.24802;
    iapp(79) = -6.54378;
    iapp(80) = -9.74914;
    iapp(81) = -1.38554;
    iapp(82) = 4.37422;
    iapp(83) = 0.769219;
    iapp(84) = 4.39207;
    iapp(85) = -3.88913;
    iapp(86) = 2.66762;
    iapp(87) = 1.93197;
    iapp(88) = 2.61223;
    iapp(89) = -3.17728;
    iapp(90) = -2.08136;
    iapp(91) = -1.50227;
    iapp(92) = -7.49413;
    iapp(93) = -6.292;
    iapp(94) = 1.53645;
    iapp(95) = -7.09815;
    iapp(96) = -0.705741;
    iapp(97) = 1.02374;
    iapp(98) = -6.51036;
    iapp(99) = -0.490585;
    iapp(100) = 4.55962;
 
  % ASCendant Sorted
  elseif iappOrder == "ASC"
    iapp(1) = -9.74914;
    iapp(2) = -9.62691;
    iapp(3) = -9.56374;
    iapp(4) = -9.28953;
    iapp(5) = -9.07758;
    iapp(6) = -8.86425;
    iapp(7) = -8.78277;
    iapp(8) = -8.71548;
    iapp(9) = -8.66329;
    iapp(10) = -8.53282;
    iapp(11) = -8.4962;
    iapp(12) = -8.15378;
    iapp(13) = -8.10526;
    iapp(14) = -8.09061;
    iapp(15) = -7.49413;
    iapp(16) = -7.38884;
    iapp(17) = -7.20252;
    iapp(18) = -7.09815;
    iapp(19) = -7.08075;
    iapp(20) = -6.86834;
    iapp(21) = -6.84317;
    iapp(22) = -6.82577;
    iapp(23) = -6.54378;
    iapp(24) = -6.51036;
    iapp(25) = -6.292;
    iapp(26) = -5.89236;
    iapp(27) = -5.86123;
    iapp(28) = -5.77334;
    iapp(29) = -5.76464;
    iapp(30) = -5.6827;
    iapp(31) = -5.67721;
    iapp(32) = -5.59847;
    iapp(33) = -5.39247;
    iapp(34) = -5.05325;
    iapp(35) = -4.68108;
    iapp(36) = -4.54604;
    iapp(37) = -4.24802;
    iapp(38) = -4.09009;
    iapp(39) = -4.08322;
    iapp(40) = -4.03607;
    iapp(41) = -3.99533;
    iapp(42) = -3.88913;
    iapp(43) = -3.84289;
    iapp(44) = -3.72341;
    iapp(45) = -3.71471;
    iapp(46) = -3.70281;
    iapp(47) = -3.68862;
    iapp(48) = -3.41212;
    iapp(49) = -3.20933;
    iapp(50) = -3.17728;
    iapp(51) = -2.69158;
    iapp(52) = -2.59407;
    iapp(53) = -2.57576;
    iapp(54) = -2.18619;
    iapp(55) = -2.08136;
    iapp(56) = -1.73528;
    iapp(57) = -1.63961;
    iapp(58) = -1.58879;
    iapp(59) = -1.57781;
    iapp(60) = -1.50227;
    iapp(61) = -1.38554;
    iapp(62) = -0.705741;
    iapp(63) = -0.667745;
    iapp(64) = -0.490585;
    iapp(65) = -0.471816;
    iapp(66) = -0.269478;
    iapp(67) = -0.258034;
    iapp(68) = -0.170141;
    iapp(69) = 0.0715659;
    iapp(70) = 0.197913;
    iapp(71) = 0.237281;
    iapp(72) = 0.655232;
    iapp(73) = 0.690023;
    iapp(74) = 0.769219;
    iapp(75) = 0.794397;
    iapp(76) = 1.01733;
    iapp(77) = 1.02374;
    iapp(78) = 1.23707;
    iapp(79) = 1.53645;
    iapp(80) = 1.66234;
    iapp(81) = 1.8244;
    iapp(82) = 1.93197;
    iapp(83) = 2.06839;
    iapp(84) = 2.26341;
    iapp(85) = 2.61223;
    iapp(86) = 2.66762;
    iapp(87) = 2.82739;
    iapp(88) = 2.83197;
    iapp(89) = 2.92398;
    iapp(90) = 2.96518;
    iapp(91) = 3.40739;
    iapp(92) = 3.7141;
    iapp(93) = 3.81252;
    iapp(94) = 3.99792;
    iapp(95) = 4.11283;
    iapp(96) = 4.12336;
    iapp(97) = 4.15403;
    iapp(98) = 4.37422;
    iapp(99) = 4.39207;
    iapp(100) = 4.55962;
  
 %   DESendant Sorted
 elseif iappOrder == "DES"
    iapp(1) = 4.55962;
    iapp(2) = 4.39207;
    iapp(3) = 4.37422;
    iapp(4) = 4.15403;
    iapp(5) = 4.12336;
    iapp(6) = 4.11283;
    iapp(7) = 3.99792;
    iapp(8) = 3.81252;
    iapp(9) = 3.7141;
    iapp(10) = 3.40739;
    iapp(11) = 2.96518;
    iapp(12) = 2.92398;
    iapp(13) = 2.83197;
    iapp(14) = 2.82739;
    iapp(15) = 2.66762;
    iapp(16) = 2.61223;
    iapp(17) = 2.26341;
    iapp(18) = 2.06839;
    iapp(19) = 1.93197;
    iapp(20) = 1.8244;
    iapp(21) = 1.66234;
    iapp(22) = 1.53645;
    iapp(23) = 1.23707;
    iapp(24) = 1.02374;
    iapp(25) = 1.01733;
    iapp(26) = 0.794397;
    iapp(27) = 0.769219;
    iapp(28) = 0.690023;
    iapp(29) = 0.655232;
    iapp(30) = 0.237281;
    iapp(31) = 0.197913;
    iapp(32) = 0.0715659;
    iapp(33) = -0.170141;
    iapp(34) = -0.258034;
    iapp(35) = -0.269478;
    iapp(36) = -0.471816;
    iapp(37) = -0.490585;
    iapp(38) = -0.667745;
    iapp(39) = -0.705741;
    iapp(40) = -1.38554;
    iapp(41) = -1.50227;
    iapp(42) = -1.57781;
    iapp(43) = -1.58879;
    iapp(44) = -1.63961;
    iapp(45) = -1.73528;
    iapp(46) = -2.08136;
    iapp(47) = -2.18619;
    iapp(48) = -2.57576;
    iapp(49) = -2.59407;
    iapp(50) = -2.69158;
    iapp(51) = -3.17728;
    iapp(52) = -3.20933;
    iapp(53) = -3.41212;
    iapp(54) = -3.68862;
    iapp(55) = -3.70281;
    iapp(56) = -3.71471;
    iapp(57) = -3.72341;
    iapp(58) = -3.84289;
    iapp(59) = -3.88913;
    iapp(60) = -3.99533;
    iapp(61) = -4.03607;
    iapp(62) = -4.08322;
    iapp(63) = -4.09009;
    iapp(64) = -4.24802;
    iapp(65) = -4.54604;
    iapp(66) = -4.68108;
    iapp(67) = -5.05325;
    iapp(68) = -5.39247;
    iapp(69) = -5.59847;
    iapp(70) = -5.67721;
    iapp(71) = -5.6827;
    iapp(72) = -5.76464;
    iapp(73) = -5.77334;
    iapp(74) = -5.86123;
    iapp(75) = -5.89236;
    iapp(76) = -6.292;
    iapp(77) = -6.51036;
    iapp(78) = -6.54378;
    iapp(79) = -6.82577;
    iapp(80) = -6.84317;
    iapp(81) = -6.86834;
    iapp(82) = -7.08075;
    iapp(83) = -7.09815;
    iapp(84) = -7.20252;
    iapp(85) = -7.38884;
    iapp(86) = -7.49413;
    iapp(87) = -8.09061;
    iapp(88) = -8.10526;
    iapp(89) = -8.15378;
    iapp(90) = -8.4962;
    iapp(91) = -8.53282;
    iapp(92) = -8.66329;
    iapp(93) = -8.71548;
    iapp(94) = -8.78277;
    iapp(95) = -8.86425;
    iapp(96) = -9.07758;
    iapp(97) = -9.28953;
    iapp(98) = -9.56374;
    iapp(99) = -9.62691;
    iapp(100) = -9.74914;
  end
 
 return;

end