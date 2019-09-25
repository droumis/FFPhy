% specifies the values for the Watson U^2 table for between 0.05 and 0.001 significance
% levels from 

watsonUalphas = [0.05 0.02 0.01 0.005 0.002 0.001];

watsonUsqr = ones(length(watsonUalphas),11,13) * 1e10;

watsonUsqr(1,5,5:12) =   [.225 .242 .200 .215 .191 .196 .190 .186];
watsonUsqr(1,6,6:12) =   [.206 .194 .196 .193 .190 .187 .183];
watsonUsqr(1,7,7:12) =   [.199 .182 .182 .187 .184 .186];
watsonUsqr(1,8,8:12) =   [.184 .186 .185 .184 .185];
watsonUsqr(1,9,9:12) =   [.187 .186 .185 .185];
watsonUsqr(1,10,10:12) = [.185 .186 .185];
watsonUsqr(1,11,1:13) = .187;

watsonUsqr(2,5,7:12) =   [.257 .269 .280 .241 .229 .226];
watsonUsqr(2,6,6:12) =   [.264 .282 .246 .232 .231 .225 .226];
watsonUsqr(2,7,7:12) =   [.251 .225 .222 .227 .221 .226];
watsonUsqr(2,8,8:12) =   [.226 .226 .222 .225 .223];
watsonUsqr(2,9,9:12) =   [.225 .226 .225 .226];
watsonUsqr(2,10,10:12) = [.225 .224 .225];
watsonUsqr(2,11,1:13) = .233;

watsonUsqr(3,5,9:12) =   [.280 .289 .297 .261];
watsonUsqr(3,6,7:12) =   [.282 .298 .262 .248 .262 .259];
watsonUsqr(3,7,7:12) =   [.304 .272 .255 .262 .253 .252];
watsonUsqr(3,8,8:12) =   [.250 .258 .249 .252 .252];
watsonUsqr(3,9,9:12) =   [.266 .254 .255 .254];
watsonUsqr(3,10,10:12) = [.255 .255 .255];
watsonUsqr(3,11,1:13) = .268;

watsonUsqr(4,5,10:12) =  [.298 .297 .304];
watsonUsqr(4,6,8:12) =   [.298 .311 .323 .289 .275];
watsonUsqr(4,7,7:12) =   [.304 .322 .291 .277 .281 .276];
watsonUsqr(4,8,8:12) =   [.296 .283 .280 .280 .281];
watsonUsqr(4,9,9:12) =   [.286 .287 .281 .280];
watsonUsqr(4,10,10:12) = [.283 .279 .282];
watsonUsqr(4,11,1:13) = .304;

watsonUsqr(5,6,10:12) =  [.323 .333 .343];
watsonUsqr(5,7,9:12) =   [.339 .353 .323 .308];
watsonUsqr(5,8,8:12) =   [.344 .363 .336 .319 .317];
watsonUsqr(5,9,9:12) =   [.340 .321 .317 .316];
watsonUsqr(5,10,10:12) = [.317 .317 .316];
watsonUsqr(5,11,1:13) = .350;

watsonUsqr(6,6,12) =     [.343];
watsonUsqr(6,7,10:12) =  [.353 .366 .377];
watsonUsqr(6,8,9:12) =   [.363 .380 .353 .340];
watsonUsqr(6,9,9:12) =   [.384 .361 .341 .340];
watsonUsqr(6,10,10:12) = [.345 .341 .341];
watsonUsqr(6,11,1:13) = .350;

save /home/loren/matlab/stats/watsonu.mat watsonUalphas watsonUsqr
