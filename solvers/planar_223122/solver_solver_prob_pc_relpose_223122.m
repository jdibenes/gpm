function [q] = solver_solver_prob_pc_relpose_223122(p1, p2)
[C0,C1] = setup_elimination_template(p1, p2);
b_ind = [61,127,193,259,325];
b = zeros(65,5);
b(b_ind) = -1;
alpha = C0' \ b;
RR = [alpha'*C1;eye(18)];
AM_ind = [19,10,9,1,11,12,2,16,15,3,17,18,4,20,21,22,23,5];
AM = RR(AM_ind,:);
[V,D] = eig(AM);
V = V ./ (ones(size(V,1),1)*V(1,:));
sols(1,:) = V(2,:);
sols(2,:) = V(8,:);
sols(3,:) = diag(D).';
q1 = sols(1, :);
q2 = sols(2, :);
q3 = sols(3, :);
nsols = numel(q1);
q = cell(1, nsols);
for isol = 1:nsols
	q{isol} = nan(3, 1);
	q{isol}(1) = q1(isol); 
	q{isol}(2) = q2(isol); 
	q{isol}(3) = q3(isol); 
end

% Action = ; q3
% Quotient ring basis (V) =  q1*q2, q1*q2*q3, q1*q3, q1*q3^2, q1*q3^3, q2, q2^2, q2^2*q3, q2*q3, q2*q3^2, q2*q3^3, q3, q3^2, q3^3, q3^4, q3^5
% Available monomials (RR*V) = q3^2, q1*q3^4, q2^2*q3^2, q2*q3^4, q3^6, 1, q1, q1*q2, q1*q2*q3, q1*q3, q1*q3^2, q1*q3^3, q2, q2^2, q2^2*q3, q2*q3, q2*q3^2, q2*q3^3, q3, q3^2, q3^3, q3^4, q3^5
function coeffs = compute_coeffs(p1, p2)
p111 = p1(1);
p121 = p1(2);
p131 = p1(3);
p112 = p1(4);
p122 = p1(5);
p132 = p1(6);
p113 = p1(7);
p123 = p1(8);
p133 = p1(9);
p211 = p2(1);
p221 = p2(2);
p231 = p2(3);
p212 = p2(4);
p222 = p2(5);
p232 = p2(6);
p213 = p2(7);
p223 = p2(8);
p233 = p2(9);
t2 = p211*p213;
t3 = p212*p213;
t4 = p211*p222;
t5 = p212*p221;
t6 = p211*p223;
t7 = p213*p221;
t8 = p212*p223;
t9 = p213*p222;
t10 = p211*p232;
t11 = p212*p231;
t12 = p211*p233;
t13 = p213*p231;
t14 = p221*p223;
t15 = p212*p233;
t16 = p213*p232;
t17 = p222*p223;
t18 = p221*p232;
t19 = p222*p231;
t20 = p221*p233;
t21 = p223*p231;
t22 = p222*p233;
t23 = p223*p232;
t24 = p231*p233;
t25 = p232*p233;
t98 = -p112;
t99 = -p122;
t100 = -p132;
t26 = t2*2.0;
t27 = t3*2.0;
t28 = t4*2.0;
t29 = t5*2.0;
t30 = t6*2.0;
t31 = t7*2.0;
t32 = t8*2.0;
t33 = t9*2.0;
t34 = t10*2.0;
t35 = t11*2.0;
t36 = t12*2.0;
t37 = t13*2.0;
t38 = t14*2.0;
t39 = t15*2.0;
t40 = t16*2.0;
t41 = t17*2.0;
t42 = t18*2.0;
t43 = t19*2.0;
t44 = t20*2.0;
t45 = t21*2.0;
t46 = t22*2.0;
t47 = t23*2.0;
t48 = t24*2.0;
t49 = t25*2.0;
t50 = p111*t6;
t51 = p111*t7;
t52 = p121*t2;
t53 = p112*t8;
t54 = p112*t9;
t55 = p122*t3;
t56 = p111*t12;
t57 = p111*t13;
t58 = p111*t14;
t59 = p121*t6;
t60 = p121*t7;
t61 = p131*t2;
t62 = p112*t15;
t63 = p112*t16;
t64 = p112*t17;
t65 = p122*t8;
t66 = p122*t9;
t67 = p132*t3;
t68 = p111*t20;
t69 = p111*t21;
t70 = p121*t12;
t71 = p121*t13;
t72 = p131*t6;
t73 = p131*t7;
t74 = p112*t22;
t75 = p112*t23;
t76 = p122*t15;
t77 = p122*t16;
t78 = p132*t8;
t79 = p132*t9;
t80 = p111*t24;
t81 = p121*t20;
t82 = p121*t21;
t83 = p131*t12;
t84 = p131*t13;
t85 = p131*t14;
t86 = p112*t25;
t87 = p122*t22;
t88 = p122*t23;
t89 = p132*t15;
t90 = p132*t16;
t91 = p132*t17;
t92 = p121*t24;
t93 = p131*t20;
t94 = p131*t21;
t95 = p122*t25;
t96 = p132*t22;
t97 = p132*t23;
t101 = -t5;
t103 = -t11;
t109 = -t19;
t151 = p111+t98;
t152 = p121+t99;
t153 = p131+t100;
t102 = -t29;
t104 = -t35;
t106 = -t38;
t108 = -t41;
t110 = -t43;
t112 = -t48;
t114 = -t49;
t115 = p111*t30;
t116 = p111*t31;
t117 = p112*t32;
t118 = p112*t33;
t119 = p111*t36;
t120 = p111*t37;
t121 = p121*t30;
t122 = p121*t31;
t123 = p112*t39;
t124 = p112*t40;
t125 = p122*t32;
t126 = p122*t33;
t127 = p111*t44;
t128 = p111*t45;
t129 = p121*t36;
t130 = p121*t37;
t131 = p131*t30;
t132 = p131*t31;
t133 = p112*t46;
t134 = p112*t47;
t135 = p122*t39;
t136 = p122*t40;
t137 = p132*t32;
t138 = p132*t33;
t139 = p121*t44;
t140 = p121*t45;
t141 = p131*t36;
t142 = p131*t37;
t143 = p122*t46;
t144 = p122*t47;
t145 = p132*t39;
t146 = p132*t40;
t147 = p131*t44;
t148 = p131*t45;
t149 = p132*t46;
t150 = p132*t47;
t154 = -t51;
t156 = -t54;
t158 = -t57;
t159 = -t60;
t162 = -t63;
t163 = -t66;
t166 = -t69;
t167 = -t71;
t168 = -t72;
t169 = -t73;
t171 = t71*-2.0;
t172 = t72*-2.0;
t173 = t73*-2.0;
t174 = -t75;
t175 = -t77;
t176 = -t78;
t177 = -t79;
t179 = t77*-2.0;
t180 = t78*-2.0;
t181 = t79*-2.0;
t182 = -t82;
t183 = -t84;
t184 = t81*-2.0;
t186 = t83*-2.0;
t188 = -t88;
t189 = -t90;
t190 = t87*-2.0;
t192 = t89*-2.0;
t194 = -t94;
t195 = t94*-2.0;
t196 = -t97;
t197 = t97*-2.0;
t198 = t6+t31;
t199 = t7+t30;
t200 = t8+t33;
t201 = t9+t32;
t202 = t12+t37;
t203 = t13+t36;
t204 = t15+t40;
t205 = t16+t39;
t206 = t20+t45;
t207 = t21+t44;
t208 = t22+t47;
t209 = t23+t46;
t210 = t4+t101;
t213 = t30+t31;
t214 = t32+t33;
t215 = t10+t103;
t218 = t36+t37;
t219 = t39+t40;
t220 = t18+t109;
t223 = t44+t45;
t224 = t46+t47;
t211 = t2+t106;
t212 = t3+t108;
t216 = t2+t112;
t217 = t3+t114;
t221 = t14+t112;
t222 = t17+t114;
t225 = t28+t102;
t227 = t26+t106;
t229 = t27+t108;
t230 = t34+t104;
t232 = t26+t112;
t234 = t27+t114;
t235 = t42+t110;
t237 = t38+t112;
t239 = t41+t114;
t240 = t50+t154;
t241 = t53+t156;
t242 = t56+t158;
t243 = t59+t159;
t244 = t62+t162;
t245 = t65+t163;
t246 = t68+t166;
t247 = t68+t169;
t248 = t69+t167;
t249 = t70+t167;
t250 = t70+t168;
t251 = t72+t169;
t252 = t74+t174;
t253 = t74+t177;
t254 = t75+t175;
t255 = t76+t175;
t256 = t76+t176;
t257 = t78+t177;
t258 = t81+t182;
t259 = t83+t183;
t260 = t87+t188;
t261 = t89+t189;
t262 = t93+t194;
t263 = t96+t196;
t274 = t147+t195;
t275 = t149+t197;
t308 = t151*t220;
t309 = t152*t215;
t310 = t153*t210;
t276 = p133*t240;
t277 = p133*t241;
t278 = p123*t242;
t279 = p123*t244;
t280 = p113*t246;
t281 = p113*t250;
t282 = p133*t243;
t283 = p113*t252;
t284 = p113*t256;
t285 = p133*t245;
t286 = p123*t247;
t287 = p123*t249;
t288 = p123*t253;
t289 = p123*t255;
t290 = p113*t258;
t291 = p133*t248;
t292 = p133*t251;
t293 = p113*t260;
t294 = p133*t254;
t295 = p133*t257;
t296 = p123*t259;
t297 = p123*t261;
t298 = p113*t262;
t299 = p113*t263;
t311 = t151*t225;
t312 = t151*t230;
t313 = t152*t225;
t314 = t152*t235;
t315 = t153*t230;
t316 = t153*t235;
t317 = -t308;
t318 = -t309;
t319 = -t310;
t300 = -t278;
t301 = -t279;
t302 = -t280;
t303 = -t283;
t304 = -t292;
t305 = -t295;
t306 = -t296;
t307 = -t297;
t320 = -t312;
coeffs = [
t276-t281+t300+t302-p133*(t50+t52)+p123*(t56+t61)
t282-p113*(t81-p131*(t14-t26))+p123*(t68+p131*t199)-p133*(p121*t199+p111*(t14-t26))-p113*(-t119+t120+t258)-p123*(t127-t128+t249)
t306+p113*(t94-p121*(t24-t26))-p133*(t69+p121*t203)+p123*(p131*t203+p111*(t24-t26))-p113*(t115-t116+t262)+p133*(-t127+t128+t251)
-p113*(t121+t141)+p123*t115+p133*t119-p123*(t115-t116)-p133*(t119-t120)
t276+p113*(t129+t171+t246)+p123*(t140+t184+t242)+p133*(p111*t198+p121*t211)-p113*(t70+p131*t198)+p123*(t56-p131*t211)
p123*(t94*2.0-t115+t116-t147)+p113*(p121*t213-p131*t218)-p123*(p111*t213-p131*t223)+p133*(p111*t218-p121*t223)-p113*(t121-t122-t141+t142)+p133*(t119-t120-t139+t140)
p133*(t127+t129)+p133*(t71*2.0-t127+t128-t129)+p113*(t115-t116)-p123*(t121-t122)-p123*(t141+p111*t227)-p113*(t147-p121*t227)
t300-p133*(t240+t274)+p113*(t132+t172+t246)-p123*(p111*t202+p131*t216)+p113*(t72+p121*t202)-p133*(t50-p121*t216)
p123*(t128+t131)-p123*(t72*2.0-t127+t128-t132)+p113*(t119-t120)-p133*(t141-t142)-p133*(t121+p111*t232)-p113*(t140-p131*t232)
-t276+t278+t281+t302+p133*(t50-t52)-p123*(t56-t61)
t282+t286+t287+t290+p133*(t58+t60)-p113*(t81+t85)
t298-p113*(p121*(t24+t106)+p131*t207)-p123*(t84-p111*(t24+t106))+p133*(t71+p111*t207)+p123*(-t121+t122+t259)+p133*(t129-t130+t251)
-p123*(t116+t147)+p113*t122+p133*t139+p113*(t121-t122)-p133*(t139-t140)
t290-p123*(t131+t173+t249)-p133*(t142+t186+t243)+p113*(p121*t206+p131*t221)-p123*(t73+p111*t206)+p133*(t60-p111*t221)
p113*(t130+t132)-p113*(t73*2.0-t129+t130-t131)+p123*(t139-t140)-p133*(t147-t148)-p133*(t116+p121*t237)-p123*(t120-p131*t237)
-t282-t286+t287-t290+p133*(t58+t159)+p113*(t81-t85)
-t291+t298+t304+t306-p123*(t80+t84)+p113*(t92+t94)
-p133*(t120+t140)+p113*t142+p123*t148+p113*(t141-t142)+p123*(t147-t148)
t291+t296-t298+t304-p123*(t80+t183)+p113*(t92+t194)
t277-t284+t301+t303-p133*(t53+t55)+p123*(t62+t67)
t285-p113*(t87+t100*(t17-t27))-p113*(t124+t260+t39*t98)-p123*(t133+t255+t47*t98)+p123*(t74+p132*t201)-p133*(p122*t201+p112*(t17-t27))
t307+p113*(t97+t99*(t25-t27))-p113*(t117+t263+t33*t98)+p133*(t134+t257+t46*t98)-p133*(t75+p122*t205)+p123*(p132*t205+p112*(t25-t27))
-p113*(t125+t145)+p123*t117+p133*t123-p123*(t117+t33*t98)-p133*(t123+t40*t98)
t277+p113*(t135+t179+t252)+p123*(t144+t190+t244)+p133*(p112*t200+p122*t212)-p113*(t76+p132*t200)+p123*(t62+t100*t212)
p123*(t97*2.0+t118+t32*t98+t46*t100)-p113*(t125+t146+t33*t99+t39*t100)+p133*(t123+t144+t40*t98+t46*t99)+p113*(p122*t214+t100*t219)-p123*(p112*t214+t100*t224)+p133*(p112*t219+t99*t224)
p133*(t77*2.0+t134+t39*t99+t46*t98)+p133*(t133+t135)-p123*(t145+p112*t229)+p113*(t117+t33*t98)-p123*(t125+t33*t99)-p113*(t149+t99*t229)
t301-p133*(t241+t275)+p113*(t138+t180+t252)-p123*(p112*t204+p132*t217)+p113*(t78+p122*t204)-p133*(t53+t99*t217)
-p123*(t78*2.0+t134+t33*t100+t46*t98)+p123*(t134+t137)-p133*(t125+p112*t234)+p113*(t123+t40*t98)-p133*(t145+t40*t100)-p113*(t144+t100*t234)
-t277+t279+t284+t303+p133*(t53-t55)-p123*(t62-t67)
t285+t288+t289+t293+p133*(t64+t66)-p113*(t87+t91)
t299+p123*(t126+t261+t32*t99)+p133*(t135+t257+t40*t99)-p113*(p122*(t25+t108)+p132*t209)-p123*(t90+t98*(t25+t108))+p133*(t77+p112*t209)
-p123*(t118+t149)+p113*t126+p133*t143+p113*(t125+t33*t99)-p133*(t143+t47*t99)
t293-p123*(t137+t181+t255)-p133*(t146+t192+t245)+p113*(p122*t208+p132*t222)-p123*(t79+p112*t208)+p133*(t66+t98*t222)
-p113*(t79*2.0+t136+t32*t100+t39*t99)+p113*(t136+t138)-p133*(t118+p122*t239)+p123*(t143+t47*t99)-p133*(t149+t47*t100)-p123*(t124+t100*t239)
-t285-t288+t289-t293+p133*(t64+t163)+p113*(t87-t91)
-t294+t299+t305+t307-p123*(t86+t90)+p113*(t95+t97)
-p133*(t124+t144)+p113*t146+p123*t150+p113*(t145+t40*t100)+p123*(t149+t47*t100)
t294+t297-t299+t305-p123*(t86+t189)+p113*(t95+t196)
t308+t309+t319
t314+t320
t311+t316
t313+t315
t317+t318+t319
t313-t315
-t311+t316
t309+t310+t317
-t314+t320
t308+t310+t318];
function [C0,C1] = setup_elimination_template(p1, p2)
[coeffs] = compute_coeffs(p1, p2);
coeffs0_ind = [1,20,1,2,20,21,1,2,5,20,21,24,2,5,11,21,24,30,39,5,11,24,30,40,11,30,43,1,20,1,2,3,20,21,22,1,2,3,5,6,20,21,22,24,25,39,2,3,5,6,...
11,12,21,22,24,25,30,31,39,40,5,6,11,12,24,25,30,31,40,41,43,11,12,30,31,43,44,1,3,20,22,39,1,2,3,6,8,20,21,22,25,27,39,40,2,3,5,6,8,12,...
14,21,22,24,25,27,31,33,39,40,41,43,5,6,8,11,12,14,24,25,27,30,31,33,40,41,43,44,11,12,14,30,31,33,43,44,46,1,3,8,20,22,27,39,41,2,3,6,8,14,...
17,21,22,25,27,33,36,39,40,41,44,5,6,8,12,14,17,24,25,27,31,33,36,40,41,43,44,46,11,12,14,17,30,31,33,36,43,44,46,3,8,17,22,27,36,39,41,46,6,8,...
14,17,25,27,33,36,40,41,44,46,12,14,17,31,33,36,43,44,46,8,17,27,36,41,46,14,17,33,36,44,46,1,20,1,2,4,20,21,23,1,2,4,5,7,20,21,23,24,26,39,...
2,4,5,7,11,13,21,23,24,26,30,32,39,40,5,7,11,13,24,26,30,32,40,42,43,11,13,30,32,43,45,1,3,4,20,22,23,39,1,2,3,4,6,7,9,20,21,22,23,25,...
26,28,39,40,2,3,4,5,6,7,9,12,13,15,21,22,23,24,25,26,28,31,32,34,39,40,41,42,43,5,6,7,9,11,12,13,15,24,25,26,28,30,31,32,34,40,41,42,43,44,...
45,11,12,13,15,30,31,32,34,43,44,45,47,1,3,4,8,9,20,22,23,27,28,39,41,42,2,3,4,6,7,8,9,14,15,18,21,22,23,25,26,27,28,33,34,37,39,40,41,42,...
44,45,5,6,7,8,9,12,13,14,15,18,24,25,26,27,28,31,32,33,34,37,40,41,42,43,44,45,46,47,11,12,13,14,15,18,30,31,32,33,34,37,43,44,45,46,47,3,4,8,...
9,17,18,22,23,27,28,36,37,39,41,42,46,47,6,7,8,9,14,15,17,18,25,26,27,28,33,34,36,37,40,41,42,44,45,46,47,12,13,14,15,17,18,31,32,33,34,36,37,43,...
44,45,46,47,1,4,20,23,39,1,2,4,7,10,20,21,23,26,29,39,40,2,4,5,7,10,13,16,21,23,24,26,29,32,35,39,40,42,43,5,7,10,11,13,16,24,26,29,30,32,...
35,40,42,43,45,11,13,16,30,32,35,43,45,48,1,3,4,9,10,20,22,23,28,29,39,41,42,2,3,4,6,7,9,10,15,16,19,21,22,23,25,26,28,29,34,35,38,39,40,41,...
42,44,45,5,6,7,9,10,12,13,15,16,19,24,25,26,28,29,31,32,34,35,38,40,41,42,43,44,45,47,48,11,12,13,15,16,19,30,31,32,34,35,38,43,44,45,47,48,3,4,...
8,9,10,18,19,22,23,27,28,29,37,38,39,41,42,46,47,48,1,4,10,20,23,29,39,42,2,4,7,10,16,21,23,26,29,35,39,40,42,45,5,7,10,13,16,24,26,29,32,35,...
40,42,43,45,48,11,13,16,30,32,35,43,45,48,3,4,9,10,19,22,23,28,29,38,39,41,42,47,48,4,10,23,29,39,42,48,6,7,8,9,10,14,15,16,18,19,25,26,27,28,...
29,33,34,35,37,38,40,41,42,44,45,46,47,48,8,9,17,18,27,28,36,37,41,42,46,47,12,13,14,15,16,18,19,31,32,33,34,35,37,38,43,44,45,46,47,48,14,15,17,18,...
33,34,36,37,44,45,46,47,17,36,46];
coeffs1_ind = [48,10,29,42,48,7,10,16,26,29,35,40,42,45,48,6,7,9,10,15,16,19,25,26,28,29,34,35,38,40,41,42,44,45,47,48,9,10,19,28,29,38,41,42,47,48,8,9,10,18,...
19,27,28,29,37,38,41,42,46,47,48,8,9,10,17,18,19,27,28,29,36,37,38,41,42,46,47,48,16,35,45,48,13,16,32,35,43,45,48,12,13,15,16,19,31,32,34,35,38,43,...
44,45,47,48,15,16,19,34,35,38,44,45,47,48,14,15,16,18,19,33,34,35,37,38,44,45,46,47,48,14,15,16,17,18,19,33,34,35,36,37,38,44,45,46,47,48,19,38,47,48,...
18,19,37,38,46,47,48,17,18,19,36,37,38,46,47,48,17,18,19,36,37,38,46,47,48,17,18,36,37,46,47];
C0_ind = [19,38,81,84,100,103,140,146,149,159,165,168,205,211,214,224,230,233,248,270,276,289,295,313,335,354,378,408,427,470,473,474,489,492,493,529,535,536,538,539,548,554,555,557,558,582,594,595,600,601,...
603,604,613,614,619,620,622,623,637,647,659,660,665,666,678,679,684,685,702,703,712,724,725,743,744,767,768,793,798,812,817,845,852,858,860,863,864,871,877,879,882,883,905,910,917,919,923,925,926,928,...
929,936,938,942,944,945,947,948,960,970,972,975,982,984,985,988,990,991,1001,1003,1004,1007,1009,1010,1025,1027,1035,1037,1047,1049,1050,1066,1068,1069,1090,1092,1093,1109,1118,1123,1128,1137,1142,1162,1170,1174,1177,1183,1185,1188,...
1189,1193,1196,1202,1204,1207,1208,1217,1227,1230,1235,1239,1242,1244,1248,1250,1251,1258,1261,1263,1267,1269,1270,1282,1285,1292,1295,1297,1304,1307,1309,1310,1323,1326,1328,1329,1347,1350,1352,1369,1378,1383,1388,1397,1402,1408,1422,1430,1434,1437,...
1443,1445,1453,1456,1462,1464,1473,1477,1487,1490,1499,1502,1504,1518,1521,1523,1538,1542,1545,1564,1573,1583,1592,1603,1617,1629,1632,1648,1651,1668,1672,1707,1726,1769,1772,1774,1788,1791,1793,1828,1834,1836,1837,1839,1847,1853,1855,1856,1858,1881,...
1893,1895,1899,1901,1902,1904,1912,1914,1918,1920,1921,1923,1936,1946,1958,1960,1964,1966,1977,1979,1983,1985,2001,2003,2011,2023,2025,2042,2044,2066,2068,2092,2097,2098,2111,2116,2117,2144,2151,2157,2159,2160,2162,2163,2164,2170,2176,2178,2179,2181,...
2182,2183,2204,2209,2216,2218,2219,2222,2224,2225,2226,2227,2228,2229,2235,2237,2238,2241,2243,2244,2245,2246,2247,2248,2259,2269,2271,2272,2274,2281,2283,2284,2285,2287,2289,2290,2291,2300,2302,2303,2304,2306,2308,2309,2310,2324,2326,2327,2334,2336,...
2337,2346,2348,2349,2350,2365,2367,2368,2369,2389,2391,2392,2393,2408,2417,2418,2422,2423,2427,2436,2437,2441,2442,2461,2469,2470,2473,2476,2477,2482,2483,2484,2485,2487,2488,2489,2492,2495,2496,2501,2502,2503,2504,2506,2507,2508,2516,2526,2529,2530,...
2534,2535,2538,2541,2542,2543,2544,2547,2548,2549,2550,2551,2557,2560,2561,2562,2563,2566,2567,2568,2569,2570,2581,2584,2585,2591,2594,2595,2596,2597,2603,2606,2607,2608,2609,2610,2622,2625,2626,2627,2628,2629,2646,2649,2650,2651,2652,2668,2669,2677,...
2678,2682,2683,2687,2688,2696,2697,2701,2702,2707,2721,2722,2729,2730,2733,2734,2736,2737,2742,2743,2744,2745,2752,2753,2755,2756,2761,2762,2763,2764,2772,2776,2777,2786,2787,2789,2790,2798,2799,2801,2802,2803,2804,2817,2818,2820,2821,2822,2823,2837,...
2841,2842,2844,2845,2871,2877,2890,2896,2923,2930,2936,2939,2942,2944,2949,2955,2958,2961,2963,2983,2988,2995,2998,3001,3004,3006,3007,3009,3014,3017,3020,3023,3025,3026,3028,3038,3048,3051,3053,3060,3063,3065,3066,3069,3071,3079,3082,3084,3085,3088,...
3090,3103,3106,3113,3116,3125,3128,3130,3144,3147,3149,3168,3171,3173,3187,3196,3197,3202,3203,3206,3215,3216,3221,3222,3240,3248,3249,3252,3255,3256,3261,3262,3264,3265,3267,3268,3269,3271,3274,3275,3280,3281,3283,3284,3286,3287,3288,3295,3305,3308,...
3309,3313,3314,3317,3320,3321,3323,3324,3326,3327,3329,3330,3331,3336,3339,3340,3342,3343,3345,3346,3348,3349,3350,3360,3363,3364,3370,3373,3374,3376,3377,3382,3385,3386,3388,3389,3390,3401,3404,3405,3407,3408,3409,3425,3428,3429,3431,3432,3447,3448,...
3456,3457,3458,3462,3463,3466,3467,3475,3476,3477,3481,3482,3486,3500,3501,3508,3509,3510,3511,3521,3527,3530,3540,3546,3564,3573,3576,3580,3586,3589,3592,3595,3599,3605,3608,3611,3619,3629,3633,3638,3641,3645,3648,3651,3654,3660,3664,3667,3670,3673,...
3684,3688,3694,3698,3701,3706,3710,3713,3725,3729,3732,3749,3753,3756,3771,3772,3781,3782,3787,3790,3791,3800,3801,3806,3810,3824,3825,3833,3834,3836,3846,3855,3865,3874,3889,3898,3902,3903,3905,3906,3907,3911,3912,3913,3914,3915,3921,3922,3924,3925,...
3926,3930,3931,3932,3933,3934,3941,3945,3946,3955,3956,3958,3959,3960,3968,3969,3977,3978,3987,3988,3996,3997,4007,4008,4021,4022,4032,4033,4035,4036,4037,4038,4039,4051,4052,4054,4055,4056,4057,4058,4071,4075,4076,4078,4079,4080,4098,4099,4101,4102,...
4117,4118,4120,4121,4137,4138,4141,4142,4164,4183,4203];
C1_ind = [39,66,85,104,119,131,135,141,150,154,160,169,174,184,188,196,197,200,201,206,207,209,215,216,219,220,225,226,228,235,239,240,249,250,253,254,261,262,271,280,281,290,299,300,314,315,326,327,328,336,...
337,345,346,347,355,356,365,366,379,380,381,392,393,394,401,402,403,411,412,413,420,421,422,431,432,445,446,447,456,475,494,499,521,525,540,544,559,564,568,586,587,590,591,593,605,606,609,610,612,625,...
629,630,633,634,651,652,655,670,671,674,689,690,694,695,716,717,718,720,721,735,736,737,739,740,755,756,759,760,761,782,783,784,785,786,787,801,802,803,804,805,806,821,822,825,826,827,846,865,884,885,...
911,912,930,931,949,950,951,976,977,978,995,996,997,1015,1016,1017,1042,1043,1044,1061,1062,1063,1081,1082,1083,1108,1109,1127,1128,1147,1148];
C0 = zeros(65,65);
C1 = zeros(65,18);
C0(C0_ind) = coeffs(coeffs0_ind);
C1(C1_ind) = coeffs(coeffs1_ind);
