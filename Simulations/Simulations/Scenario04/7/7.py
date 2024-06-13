#!/usr/bin/env python
from numpy import *
data_1 = [
array([29.33740229168414, 34.218847308079255, 30.83479448038041, 28.171744148592296, 33.97135455591386, 31.584251141869142, 30.65199753735197, 32.671047000538586, 28.344140497157888, 31.65066818692436, 31.696979064407653, 28.77162043945803, 26.061027148062077, 26.463885443508058, 26.095361672665767, 27.047198783190563, 26.88554383073735, 27.64544459862953, 27.365182947065684, 26.392941546723385, 25.855921799087536, 27.744079287895026, 27.981494145104307, 27.94216074256008, 27.080996060956597, 26.245149211078424, 27.573787188107925]),
array([34.0159281608529]),
array([29.850213091991918, 32.446952363234985, 31.569118852025348, 31.185518841787253, 31.22727897382604, 32.761675960829905, 31.786953982846054, 31.647650636810614]),
array([25.680940398026415, 25.486910928826198, 23.852848841748756, 27.268351640594247, 28.02537666717055, 25.53125898610899, 17.97043762863712, 16.854307823732967, 17.244426132962705, 14.226741023862115, 14.122025160164426, 12.816284131654992, 10.127663825570279, 4.750063418954519]),
array([28.226203538631843, 28.33908212006166, 29.2970759335138, 27.183465534706052, 27.346213330289494, 26.107628288562015, 26.336654955838252, 27.883442158985442, 26.685893324165974, 26.041560882694665, 25.753014409313966, 28.099522908902756, 26.016520745390626, 27.98842034663568, 26.340128902540716]),
array([25.068103984491618, 23.27327739834322, 24.896373002342205, 27.74540894509372, 24.701082225438803, 26.807553487591232, 21.550692393738682, 21.345995844442804]),
array([23.882043952778858, 27.940487204153243, 26.45663137347825, 24.67047005301128, 27.310931563290286, 24.82340091307359, 24.56507023723958, 24.37974723863319, 23.662050575794037, 26.005247791616465, 27.37761132439941, 26.780223582916893, 26.12429501468575, 27.823291719204473, 22.741129456870382, 22.9479489502373]),
array([28.243925757335404, 27.674870702150393, 28.062193611335527, 28.065374745662925]),
array([28.80214878254852, 29.59757199718006]),
array([29.627647347673687, 25.106676485518026, 25.5119507305259, 27.285683839728247, 27.003117025301776]),
array([24.143844637952157, 25.13386212428649, 26.831044618355655, 25.949198051447542, 24.00186957675205]),
array([28.47714795497584, 28.099592101947092, 27.514487618791858]),
array([25.230716982693618, 27.17715821067323, 24.043022474794185, 27.71369039670093, 26.021299597594776, 25.194219425230326, 26.846724654210373, 23.61186961559158, 24.566747224416908, 27.279531335793664, 27.326320526235214, 27.785998688606483, 25.20250253345026, 27.150556258629177, 25.347077832195293, 23.945274535915296, 23.64253194470438, 27.510103793419958, 26.415556248771452, 25.983905112736327, 26.156073393857028, 21.538098450592738, 20.909348776822398, 20.966775875617863, 20.005357294688775, 19.629331885375137, 18.852728671101165, 19.505439019733668, 19.98093827171575, 18.528326675380498, 19.523885787277166, 19.810197067857146]),
array([25.85449126677263, 25.654460840622257, 24.213165563993126, 26.405193264266064, 26.220957766275745, 25.99341058773127, 22.09848088092389, 21.53569660664482, 21.168572491605328, 18.154675530270943, 19.094894936898022, 19.231384594026313, 17.482952905858408, 16.885307110731677, 18.61252213722588, 19.720597459063175, 18.31319461613066, 18.644844682219237, 18.492370492767844, 19.367754297762357, 14.852029741495597, 15.751889840674957, 15.498396001448246, 15.917445244003524, 15.613674696748864, 15.838581574716258, 14.46498772949649, 14.46292287024443, 15.01687159905494, 13.139441270113712, 12.498936421117316, 12.931884884102514, 11.867400665253239, 9.977433039712224, 8.294030625564083, 10.785421256340932, 11.282727346667135, 10.825984853613821, 8.506285214081451]),
array([27.464147074367602, 27.128738596798854]),
array([27.267906511255475]),
array([24.597349901656223, 24.87682316793048, 26.620969283525387, 26.58517378203593, 24.992301755082824, 25.648333481685405, 27.269878418192107, 25.441844994672906, 26.135006908572297, 24.983619392330688, 25.0968012482739, 24.579132374801397, 25.23710638477012]),
array([23.293674881865392, 25.197495899456893, 25.033691459532495, 25.880621438698956, 20.76844206805628, 21.756685328152344, 16.999899760004677, 19.536042580885773, 20.26725003351108, 14.075923961409337, 15.655481408264365, 15.55288833488341, 14.655476100717141, 14.834639045305702, 13.274914574752618]),
array([23.262861808541135]),
array([24.357878387006146, 23.73442836591303, 23.80056324015976, 23.322121926055786, 23.990856537735148, 24.3760164647617, 23.998562757783006, 24.524561508451857, 24.126692393297958, 23.610895741817355, 24.510285263812182, 23.046176924469084, 21.375304783203937, 22.429492416997327, 22.615613453045935, 22.474106998133017, 21.65073957114743, 21.35708461884031, 22.13946566278357, 21.714715483921854, 22.767642004192652, 22.07095865096599, 22.208150239214113, 22.21426794865348]),
array([23.210886255315963]),
array([22.31561682966509, 19.48784546148858]),
array([17.121019905203966, 19.372745678845657, 19.664875536398398, 17.84502518105416, 16.329251958900798, 19.304694661463014, 18.026181243000686, 19.25822935762069, 17.124834418053766, 15.707299756470238, 14.020272573243567, 14.633041319374668, 14.662903051413693, 15.167809872481227, 14.666654711376683, 13.609928595109855, 13.795298648409135]),
array([20.990739675195638, 19.006014080847667, 19.276604008523584, 18.42564843226892, 18.808542009977213, 19.82247149283231, 19.16222841553789, 17.503564930697188]),
array([19.2142631224754]),
array([18.106964312146367, 17.457170958282, 19.726413160984983, 18.54021302620036, 16.430704301668506, 19.526032901748277, 16.049184554893774, 15.816257177640066, 13.446991546592773, 13.613029542904657, 11.8255289684535, 11.6017340272812]),
array([16.59730121980596, 17.266290954211005, 18.55693225512882, 19.719719712291738, 18.882955863794823, 15.322404293139828, 12.330994353957427, 12.54755566340421, 9.199761486132855, 9.41969814738523, 9.989610377993223, 10.107343810781336, 3.8184854582923524]),
array([17.58061198666617, 18.666487360238587, 17.37018748568649, 19.366291980240597, 17.299263157691037]),
array([16.532351434968845, 17.08676740222875, 19.169863175899593, 19.09397640726858, 19.731536523718702, 19.78376858629632, 17.018940078621426, 14.608238251057635, 14.682767054617688, 14.839049390516967, 15.855844944322003, 13.574364490803983, 13.659285950479596]),
array([17.98867893662437, 16.916800462352626, 18.817497697926346, 18.569290513854043, 16.388026722855766, 15.980387899560382, 19.536028908359512, 16.846369194787066, 17.73209014597162, 19.921941308992167, 17.492295923116004, 18.82868509564732, 18.54826462565635, 16.09196281730351, 16.991568683655093, 19.72354869075535, 17.236526085445142, 19.888407236941553, 16.421810390261296, 15.378203671232587, 14.158490506307496, 14.234528979353234, 14.460776530342471, 15.852678300921465, 15.83479123833987, 15.578692066711714, 14.446927321930966]),
array([17.33321309401799, 17.217823971754033, 19.294403214077956, 19.0668037591567, 17.571603168404135, 18.05973702406122, 16.584161290923642, 18.390174401646185, 17.99912999039198, 15.794625126572562, 15.776238778332937, 15.746057570901538]),
array([18.724403946588666, 17.526013989328966, 16.110782480208023, 18.267742786480873, 16.624696994223267, 18.58681736285674, 16.327123798664275, 17.49545813021288, 18.706220443702733, 16.65015952895013, 16.32185359326165, 18.64595854230487, 18.73248816942374, 18.313154104158468, 17.488211371953348, 16.12615461266533, 16.406324568786513, 18.08392069698875, 18.983376256081343, 17.235434391216423, 17.30099470905085, 18.343986867112747, 17.210855132616604, 16.02760982909448, 18.648068841366534, 18.114668571160497, 16.44596798565223, 18.601782598629484, 18.696953668725254, 18.561571585211915, 16.690846667026843, 18.701144001431846, 18.814780050835935, 17.137624190008108, 16.699821544882212, 16.732635659529187, 14.430246748867603, 14.862081594292631, 14.16289217265384, 15.810312267161658, 15.861736365212472, 14.192848111220508, 15.741170349553851, 14.468628300793092, 14.700405767923467, 15.534776621351593, 15.423552628953175, 14.416532303916135, 15.924444001780696, 13.96309991809264, 15.875643881461865, 15.63050546303023, 14.62610164161934, 15.729802181551369, 13.874881667900413, 15.060703670429909, 15.82519898135349, 15.623456178925187, 14.831255440334687, 14.142380251186257, 15.43020136397605, 15.131303453350416, 14.775167631717451, 14.603549926993075, 14.451125561918726, 15.288016253070436, 15.19166023749709]),
array([17.715241954390823]),
array([16.1434713790833, 17.736179399145794, 17.131635170360138, 17.28188271708426, 16.913933067000233, 17.64801449534616, 17.743984385481617, 16.004614663300174, 14.091753926583836, 14.021919344136396, 15.072590708935024, 13.827417997385876, 15.593372148030543, 14.199923369362502, 14.015814240218635, 14.755974031167028, 14.482367629185742, 15.8407897730681, 14.158332469819939, 12.916177333600276, 13.617082444848558, 13.451884459399977, 13.029690647717537, 12.749997499280504, 13.458898526786045]),
array([17.642043916776366, 17.460274260649353, 17.077033116252824, 17.641560285364488]),
array([17.002448486921296, 14.501199144575056, 14.726391940603863, 14.720351565287697, 15.62614046055327, 14.148302597102253, 14.832381699559159, 14.66298022329594, 12.826450601716553, 13.202244027769867, 9.268465289633966, 11.363281001584292, 10.21455797759655, 9.804235440655926, 11.009345115340395, 9.330050620410782, 1.6081057874293108, 1.140429040612458, 0.0]),
array([16.50731919482335, 17.11628327861288, 15.604323844103718, 14.619524437769844, 15.457537850770208, 14.330898571608255, 15.910619223490235, 13.186440047153546]),
array([17.073583839111233, 16.38498279916974, 16.68090411920897, 15.70249331407948, 15.516059860144729, 15.107961688284309, 14.098031816426689, 15.95537413221856, 14.774732258828418, 15.892587634221975, 15.153820020905394, 15.328947891522962, 15.617096592930487, 12.56103910569853, 10.273490836397567, 11.342430852576816, 11.068446392207562, 10.192708791852262, 8.728535462711035, 8.121735090177086, 6.185187803698185, 6.786799470128294, 6.543724615510684, 7.11130444969556]),
array([16.982548651095957, 16.75469848761725, 16.30719793040732, 14.939431493854695, 15.235084879002219, 15.406990075873837, 14.654876001919858, 15.444651236189568, 13.880177051988717, 15.433598946697042, 13.844616546209016, 15.623049894438601, 13.927707060883252, 15.83512012368033, 15.954036524670482, 15.760408837558108]),
array([16.2344326932782, 14.449570806815572, 15.469802250957473, 15.377398696749, 14.522356742766643, 15.906575663260057, 14.808434317302558, 15.109271331129726, 12.413890324020523, 12.326504492991345]),
array([16.415143473018098, 16.671312262751645, 15.828004717785264, 15.16767260710529, 15.652868369032438, 15.713621230514539, 15.358836348642447]),
array([16.09519308811405, 16.587999447324865, 16.159605871112685, 16.630485297662965, 16.020814555675436, 16.146770790022792, 16.359032986393068, 16.411346389635035, 15.234292100377685, 15.77389176347626, 15.860637293601622, 15.144373039969024, 15.602401279502823, 15.509907718343499, 15.577031588442281, 15.257960982921986, 15.60112006758173, 15.318926593588975, 15.520429283569229, 15.901467898576689, 15.531975781070285]),
array([14.130462670988164, 14.828098532349136, 15.125499328294778, 15.914870630075347, 14.462216081503176, 10.832599606207367, 11.553391379693686]),
array([15.125909618063194, 12.580932531219313, 13.767063026805282, 3.7384355950756447, 0.39293332937476205, 0.0]),
array([15.254512725552145, 15.156086903750268, 15.054471435128713, 15.228184779390464, 15.321625257943241, 14.963741850472093, 15.841152631680567, 15.461387476398292, 15.011162259134835, 14.064280297130203, 15.03236225073416, 14.761507240277021, 14.045583358504278, 14.231653404498095, 14.429931870136173, 14.626750955736572, 14.67673016762159, 13.214734199261969, 12.38399443188631, 11.7946669796945, 13.500082356044462, 12.005474971189782, 11.893957774881597]),
array([15.243745886707934, 15.312428387593755, 15.410958047289826, 15.36366652582931, 14.71742285146281, 14.725484985289114, 15.433766011380428, 14.440467794371724, 14.932504535963686, 14.812335478503664]),
array([15.396034865709895, 14.7131193615378, 15.390341750682936, 15.276428614914082]),
array([14.493894426002724, 14.45296209325232, 14.408138808736146, 13.900586870323338, 15.035933949584868, 14.227960033861361, 14.538473863200805, 14.056579106403609, 14.238935091267063, 14.056114648751613, 13.623364457110831, 11.838454215945246, 12.516139748892833, 11.96689870062811, 12.344645439527396, 12.333650395043183, 11.904136267661137, 12.149099343735413, 12.924784405639496, 12.982391447869265, 11.698546342102183, 13.241942901769537, 8.968969748379148, 7.661074405461708, 11.40410760919019, 10.105669683656593, 11.556098534405649, 9.099130604384525, 11.172807030301165, 7.504933634283437, 9.844024995222613, 8.905969891636852, 10.188979217985823, 10.157299634033182, 11.579412032765694, 7.830608075405983]),
array([13.989225092487379, 10.447530217974693, 10.541192987307994, 1.3860022387436788, 0.0]),
array([14.092301793423408, 14.896003262922333, 14.39384352151988, 14.101112775279521, 13.032602340938952, 13.519149107512192, 13.111259205748926, 12.599039692812482, 13.414874058438102, 12.929955768557981, 12.547776128681006]),
array([14.077179619813991, 13.903150058089581]),
array([14.456661100719163, 14.900227872845681, 14.117310326197426, 14.606801721072564, 14.195388300997513, 14.856606723562404, 14.522642387718742, 14.219858139628148, 14.772603933491332, 14.564323378854484, 14.75190583915691, 13.938732497714014, 14.065998946420299, 13.199614139570311, 12.604173117717107]),
array([14.20318266709107, 12.691325702117606, 10.629372662613802]),
array([14.165645700441456]),
array([13.835452587675855, 14.33157485839102, 14.023193607828812, 14.450491653357977, 14.081081513176924, 13.351864883172185, 12.776318350136473]),
array([14.133822323919956, 11.002156950006327]),
array([14.077552440065356, 13.843605808136157, 14.007368734217902, 13.365733322354773]),
array([12.830922478424702, 13.467180020423918, 10.99910150597929]),
array([13.820776080806354, 11.658703182614865, 13.439090389010982, 12.97769644890902, 10.397595402774046, 7.355790110914807, 7.7884079450142085, 5.390660138040117]),
array([12.537610760607771, 11.776108989931029, 13.159861836142953, 13.569026912527818, 12.272281120296679, 10.391956553056644, 7.951003725463778, 11.462835799838636, 9.137240770537304, 7.84437031265261, 10.679132822191608, 10.924403333283454, 8.148570573075878, 8.765932418487285, 7.359342735778066, 10.933886258337871, 10.979264934301304, 11.387189550330005, 10.208328167745346, 11.420763578293371, 8.938404523100662, 11.253559280655123, 9.63861457195442, 9.371938105582519]),
array([11.737408528331134, 12.830209488015416, 11.634363259214865, 13.114942220871045, 11.891520807687233, 13.816831413059361, 11.251434071011285, 8.596732897367481, 7.287423162299592, 9.923395570751651, 9.364314275295492, 9.984873430557675, 6.897706647043876, 6.352702079273675, 7.076894970782996, 3.7774263068828757, 4.923337382171933, 4.23187974773526, 3.4044842637935657, 3.2888253167854673, 2.6923485547629875, 3.2277454272409125, 2.3927294624356783, 2.059107443287507, 1.854692487055439, 1.6827514762499778, 0.8968822422501572, 1.066607977419542, 1.527855676697901, 1.595337126845146, 0.4889622102751915, 0.3563551558074465, 0.0]),
array([12.174723238330458, 11.825172242060733, 12.854115681322867, 11.88883559089438, 12.31093355284471, 12.785670965886297, 13.256836450997103, 12.205224709794473, 12.691871631093811, 12.405726324242577, 10.032925981608917, 10.243707237735848]),
array([12.864114366617954, 6.944898818959658, 4.001800048789406, 2.60971482994756, 2.426233117588785, 1.4357510919386494, 0.7229683028778117, 0.30991718726438455, 0.0]),
array([12.934261021891851]),
array([12.668613854339911, 12.042943972698048, 7.383215315074374, 9.93597437585179, 10.24320233327429, 6.985864093730825, 1.4784346990060855]),
array([11.596107442091077, 11.468587563302586]),
array([9.908617131977074]),
array([9.964602357857219, 10.35865810230218, 11.180165549564444]),
array([12.34549843019966]),
array([12.17089335282597, 11.944212095280987, 8.093636499236979, 11.481637528687786]),
array([12.087991270087986, 11.962170246846892, 11.208421898386659, 11.240681016814676]),
array([8.478439871147764, 8.424610277868577, 8.72291294864227, 9.24457816002749]),
array([9.728857590716954, 8.114671223951788, 11.047906120090197, 10.391569113658917, 8.202612737028408, 7.7322705134173395, 8.727180390460845, 11.411394978054334, 6.5072080093633575, 4.835301522720541]),
array([11.09019132122505]),
array([9.15075829162436, 10.451451549899298, 10.630660423737922, 10.72808874179964, 10.077202803203074, 11.19551043218271, 11.285374163478899, 9.093446798521754, 10.34734494532405, 10.312189845177576]),
array([10.848430201812457, 11.505640030456435]),
array([11.216161285798105, 11.239809076300407, 4.504629332154106, 2.071616453759843, 0.0]),
array([8.785609321392961, 8.502836971836302, 7.869954131257458, 7.761430569524306, 10.822100554323335, 11.324264612055982, 8.854603659414563, 10.30219871801659, 7.975739201017097, 9.580560979086767, 9.76293694630557, 9.82807244839346, 10.247258978152683, 9.298919983745762, 10.698020015631187, 7.774266047818061, 8.403804954423391, 8.23111558926767, 7.889644074746656, 8.860731231208677, 8.052172812374982, 8.655006887946996, 9.630219848705911, 11.302639343460207]),
array([9.742241485975756, 9.674826351051985, 9.125497357327376, 8.755874685716755, 9.425949548827726, 10.936796528647397]),
array([10.761600269053172, 10.519983379740541, 10.516907653600528]),
array([10.059044853401545, 10.607342871198156]),
array([10.338078242797145, 10.31391534170584, 3.848658492071187]),
array([7.509550500034745]),
array([9.578258934184046]),
array([0.867101083597215, 0.0]),
array([8.696048893313097, 7.735343847764165, 9.291843617541225, 8.509742283824863, 8.386362725400836, 8.82768011107796, 9.624079288128653, 8.804016008915122, 9.64703661172538, 9.091719055495442, 8.867466306537318, 9.435315408346161, 7.253110154154158, 8.084572237297527]),
array([3.4298645591775765, 2.190327750823574, 1.2749442871142902, 0.6743204432790662, 0.0]),
array([7.421142166612658, 7.419599639936845]),
array([9.045790226971574, 9.404022892910374, 7.266873942180833, 8.482284780452803, 8.061348162877618, 7.583038827140478, 8.912094767886373]),
array([8.267024926703547, 7.077730420998782, 4.113399196036071, 4.0736517061915585, 5.191962500419554]),
array([8.879213583972444, 9.072546254433956]),
array([9.461025411601584, 6.351575139092534, 5.238656705563998, 3.0575833883391965, 2.0625883322601197, 2.4725765653880982, 1.7966157921195203, 1.7861266041853707, 0.6039050461289507, 0.6633516420777953, 0.0]),
array([8.984018076294594, 8.953008955751303, 7.820843449775546, 8.346037976508729, 6.394774408312863, 6.958711173934608]),
array([8.223908649593678, 8.983480065188342, 9.097802440454439, 8.528360596530783, 9.007025877205214, 8.254070396236974]),
array([8.589645448215778, 8.809740655370627, 7.782070444789669, 5.7751667060537795]),
array([8.603974925623875, 7.706592099109178]),
array([8.182034936466255, 8.649450086749912, 7.817074350029679, 7.845432729759023, 8.626679578077644, 8.074269195589114, 8.070128391231693, 7.681841603383521, 7.100215620664954, 4.137257931713401, 2.9489565359360252, 2.90018465462744, 2.7784491935579796, 2.398378942635832, 2.3088559964183664, 1.9251098789130916, 2.280222172674349, 2.316807429976148, 2.0015376484171172, 1.5504787003812823, 1.001854296160074, 1.6287198538288385, 1.2128028781409306, 1.7704208669864616, 1.4713567459010908, 1.1271770148384954, 1.1062958863036225, 1.7817858244116025, 0.30912690143549154, 0.754962032646794, 0.14101343952206435, 0.5363392923541336, 0.0]),
array([6.144549331696681]),
array([8.07666693311793, 4.339190308508108, 0.0]),
array([5.6940342540657785, 6.0693864355998866, 5.034157577754389, 5.3223966099602125, 4.922820425512074, 2.8286192053929464, 2.7656999170922547, 2.4717934408404214, 2.3304088087012937, 2.1503625199038963, 1.6208445026906095, 1.564896306060204, 1.5389562226258007, 0.7528480221732727, 0.7415202058887779]),
array([7.238321017504388]),
array([7.445887455393516, 7.249479247634131, 7.692880895580399, 4.252366585887646, 4.703049254827936]),
array([7.102406105829404, 7.157882607321076, 4.20475250549386, 5.065057307821459, 4.208155367861503, 5.136777230714799, 3.780385515749111, 4.886791141248491, 3.5397190521269333, 3.5175194720164202, 3.464234494222702, 3.398033822519415, 2.867459259154579, 3.0513646642326258, 2.8789300873908465, 2.9888446540282154, 3.0291157007162544, 2.309832468644907, 2.1136154839081347, 2.1108225676475225, 2.4697871973810295, 2.4759288880110484, 1.9682654124294685, 2.254455257394707, 2.46366743591191, 2.5204023095715327, 2.2532490719989964, 1.4345110581883804, 0.999187894532308, 1.6944537351349827, 1.5651558174381641, 1.3061730794220152, 1.005161678072846, 1.014620179958695, 1.6891343902472191, 1.2549159001355998, 1.1384502577165376, 0.4563161926746548, 0.6302710958909479, 0.28203946390457885, 0.0]),
array([6.998585424835253, 5.322126913745456, 3.410594619750566, 2.722575158914762, 1.9231325933951124, 1.3558830440951302, 0.0]),
array([7.351980560486321, 6.237721144247056, 3.944310011287266, 4.774954491556225, 4.6221516659090645, 5.080272642308018, 2.8117096155683985, 2.8288722601690877, 2.420010655410645, 1.8038804018150958, 2.563721431269518, 2.1170370746187435, 2.034682586559422, 2.164827413544016, 2.3778506438778844, 0.9312597847949898, 1.6878406161117818, 1.2288114415180846, 1.07010833424529, 0.8003002100431335, 1.7899843260868327, 0.2862643167176358, 0.3689379129441038, 0.18751912011344263, 0.0]),
array([5.408045736333071, 6.539396696588792, 4.978247283180438, 4.639835760821214, 2.64960095123965, 2.333171248396179, 2.1479811638246744, 1.8816297073545953, 1.0626136064819751, 1.4526544740323533, 1.0406320217484084, 1.2631233161579758, 1.0051305615504011, 0.3079496562501522, 0.0]),
array([7.307131999874721, 5.73513797206917, 6.728216286524972, 5.382394125423734, 6.51267446174767, 4.826061741681103, 4.455710479334105, 4.86783568071907, 4.573785947675677, 5.279635507180707, 4.459217404927045, 4.290366106232773, 3.1966383675945487, 3.1696724797552602, 3.119484482672619, 2.2526073179523167, 2.3110306500805358, 2.2677068334431825, 1.8005230328453323, 2.3721419484338533, 2.4463242638188873, 1.6321335520399287, 0.8679506074178963, 1.4048940528241518, 0.9812407519417922, 1.5017128199583982, 1.6623851967295875, 1.450452743511304, 1.1124533309267237, 0.8626462487625001, 1.2935947816977436, 1.4599456083962785, 0.8504382531982869, 0.6759406984197286, 0.6078031611150273, 0.4688824266726643, 0.5783766587448218, 0.12565088564153407, 0.0]),
array([3.7757457064808246, 3.330719850678215]),
array([6.786811681812381, 6.793683720097787, 6.933222042225431, 4.09200966565602, 5.193825596940142, 3.078562817134933, 2.961252395882955, 3.0008333602794925, 1.0454811165660303, 1.3392141375591213, 1.1348529423329006, 1.5404470475160732, 1.2926446218638232, 1.3921841577780656, 1.549632912884452, 0.0]),
array([1.8835639935045239, 0.0]),
array([2.482134892389956, 0.0]),
array([6.757993185356417, 5.802895807813242]),
array([3.420040220739093, 1.9252393317162082, 1.7098386874080715, 0.0]),
array([5.709164144853028, 2.8602902729572794, 1.567447154078081, 0.9858323133823542, 1.025705897988812, 1.353493278447952, 0.0]),
array([6.492044772618557, 5.871151067190798, 4.775615510272292, 4.175960194178153, 4.2383177635579425]),
array([5.860661687175963, 5.237441059630406, 3.5519933696698005, 3.4590932944141635, 3.3302560474179463]),
array([5.281186883543882, 4.4722001320695375, 4.555777123887003, 5.092541406828336, 2.7566736953334554, 2.092301174854824, 1.832789191296821, 2.3321353031978758, 0.9023019615284321, 1.6389522108459227, 1.05175162900558, 1.0226813571928095, 1.731597027517723, 0.9177909456394115, 1.0111112006457463, 0.5214580620189797, 0.0]),
array([3.3429096704331243, 1.326710349714565, 0.0]),
array([5.803224120042801, 4.8718503772536845, 4.725279682030465, 3.9681482305333913, 4.881324455698565, 2.7559856677824377, 3.007797448062035, 2.8869005549357922, 1.9125855007534733, 2.511756511841349, 1.9697744024556767, 1.2961477250016673, 1.204714466869441, 1.312258810813439, 1.744330973907454, 1.3606461665351388, 0.9186466372331351, 1.2554087573110861, 0.0]),
array([3.2977306723551765, 2.324805859163378, 2.2207489531378135, 1.550684868244649, 0.874780316983406, 1.3654155034667361, 1.6124886520822428, 1.193209094525614, 1.6225270985337288, 1.164229673474332, 1.1097798313277003, 0.0]),
array([3.936186548414562, 3.876486153247739, 2.609518912085066, 2.936020390914893, 3.5179244171984885, 2.777524727493426, 1.958294424143265, 2.0847756807979083, 2.222568567330116, 2.3517715624816367, 1.60383596235849, 1.4794740471174717, 1.2008455641770355, 1.5758231715818145, 0.9627834647489132, 0.47017036580860216, 0.26764145708504405, 0.0]),
array([4.9029595970358235, 0.0]),
array([3.854643157221236, 5.190287720025494, 4.66832083504573, 5.157892483082108, 3.1233711165336815, 2.7485200258994773, 1.054743601629022, 1.6212956783012042, 1.1997104636916034, 0.0]),
array([4.737075286784357, 4.804291595566376, 4.001034214050078, 2.885403798770745, 2.6144541683360316, 1.9907693905674226, 2.2166556497144088, 1.2164357089761837, 1.4194790881183783, 1.1304856381511459, 0.42481392576501814, 0.23228813514910085, 0.0]),
array([3.876758453962805, 4.699019176928181, 2.6726363994297015, 2.9201921731789007, 2.8192639936100843, 3.5387460529425017, 2.789213656511781, 2.850033391988789, 2.389174570605503, 2.2741609683220037, 1.9591528271708936, 2.534273386068956, 2.328654346083318, 2.012635176616676, 2.102120234135672, 2.2372139814445307, 1.3622546538144533, 1.3472799031032423, 1.3983341174841033, 1.2151161471305882, 1.7988796033649108, 0.0]),
array([5.053975315718247]),
array([2.637837373481294, 1.1876459386762486, 1.6780494152864225, 0.38275446455437867, 0.0]),
array([2.4585411935894066, 2.542876107593588, 1.568457084930581, 0.9713377794001689, 0.0]),
array([3.6363749685207565, 4.462715271743531, 4.142684117923984, 2.8298185732142267, 2.9707200075969706, 2.267164472331938, 2.043127471853464, 2.37106602174586, 1.8763942372640394, 1.5685637924971934, 1.1677670742009134, 1.5419888681243679, 0.0]),
array([3.9958253155923322, 2.8060139528322865, 2.9855668195900384, 3.1076502917265967, 2.6005971642623433, 3.0958789321299087, 3.4969754361887864, 2.3383929217763466]),
array([3.5534649141608887, 3.223724178895047, 3.1345767955004087, 0.0]),
array([0.7855369967905674, 0.7760986422006538, 0.0]),
array([2.2249130565390955, 1.9913282432050334, 2.3462525760540753, 1.4857918944491402, 1.6937749002805182, 1.5082890975416414]),
array([2.9749007194541375]),
array([3.380109951145142, 2.9666619457395367, 1.8827093165759454, 2.2933592875803104, 2.5758904167842185, 0.8720241884580122, 1.5919999022357938, 1.4172489585058958, 1.608029196487431, 0.13405941640445795, 0.0]),
array([2.615247514504017, 1.8167324539677567, 1.5670487624231244, 1.2635793472440757, 1.288255733822316, 1.339672814673245, 1.597379455890405, 1.256417498457007, 0.684829125036917, 0.08151194174482593, 0.0]),
array([2.193055716763033, 0.5218861963990432, 0.7683115074928126, 0.0]),
array([1.875197718156748, 0.8354932441176821, 1.4333061484002587]),
array([2.494954827260799, 1.099428743210711, 0.43714136300797474, 0.0]),
array([2.7413968614182957, 2.424073702827854, 2.2026994848113732, 1.5650340234867968, 0.0]),
array([1.8815754801532405, 2.5082172822390305, 1.0332817095142361, 1.2129873614786888, 0.0]),
array([2.7223682745722937, 2.812693596589318, 2.159978782296417, 1.9886740336432744, 2.2754798084541434, 2.2326486499263916, 1.8516792458214657, 2.1611134752630656, 0.8649536739447825, 1.7205413907122353, 1.1309624995720373, 1.3705529977059092, 1.423494765440358, 1.460562382606314, 0.8804183859412654, 1.542602338561902, 1.1564454850400696, 1.773276151203651, 0.8183126660232742, 0.9827529179697914, 0.33709435479382654, 0.5924608414721991, 0.0]),
array([2.6262586390660516, 2.0417838837789604, 2.0550888715144464, 1.88014524071943, 2.3136617169045413, 2.1532018099470807, 2.413394103999813, 1.0973815559098345, 1.2713172471312335, 1.4331415875562064, 0.5758138185293018]),
array([2.916366535015803, 2.6623406010508592, 2.320411080906273, 1.3560340861379203, 0.0]),
array([2.6688862919428082, 2.7727243921824125, 2.737934209492335, 1.9389567303288437, 1.828297881106898, 2.3242167483827654, 1.9858707040530705, 1.1282412899782883, 1.2815976490348377, 0.889461187426054, 1.469022003023529, 0.9204027695026562, 1.4933907417012044, 1.6522367350270641, 1.6703416058295364, 1.5292473223677714, 0.8654906801782629, 0.6827790400524428, 0.1596836270469466, 0.2660292425473888, 0.6813171652014837, 0.7135496658225691, 0.7376748421680722, 0.40242259404668673, 0.0]),
array([2.8259222067230803, 2.1606556183344616, 1.5438634328952872]),
array([2.039778945166428, 1.050669017369604, 1.4876878792445356, 0.0]),
array([2.294219012249309, 1.5793342550157652, 1.4314270024106395, 1.119220432209581, 1.634610647339148, 1.2615895915239141, 1.424310147842423, 1.0258158929508157, 1.1437509429962516, 0.6166820661060952, 0.0]),
array([2.304358072226014, 1.6479125484609098, 1.4851692064925874, 0.9767251812107299, 0.8960085896780134, 1.2792073858447313, 0.9571568056876442, 1.585391726710719, 1.3862320842514506, 1.1443170004063306, 1.7717261231547317, 1.2603740888698294, 0.931050351988964, 1.585953561580371, 1.7879958689791906, 0.14488198524255758, 0.003967349477572157, 0.04378546644502927, 0.0]),
array([2.4191340746892043, 2.188867311647174]),
array([2.3098744559552395, 1.1400255218022133, 1.1268474762640732, 1.55272764574716, 1.4673795093026871, 1.44945533945317, 1.7150860302220823, 1.0417180649336863, 1.5815674904123878, 1.0134516624897996, 1.1362408789335332, 1.6943139266297904, 1.0378261582536772, 1.6534846421401195, 1.1292245940868804, 1.7209223040849626, 1.2082972625758286, 0.7266615396603981, 0.40056893897205154, 0.6041092385421796, 0.7741990536740398]),
array([1.869839081799111, 0.0]),
array([2.162872895291539, 1.9757967641965246, 2.0371689423959567, 1.3473320526028763, 1.2658919237409751, 1.0371600687240523, 1.047681321742412, 1.7536402938941313, 1.6336065930868187, 0.29889174176209843, 0.0]),
array([1.995774367736803, 0.7993289826305816, 1.3523313615541284, 0.962721504994738, 1.2388305985494301, 1.0120554706841036, 0.0]),
array([1.7278603333211713, 0.0]),
array([1.288890322761137]),
array([1.8039320697637569, 1.6269558237857626]),
array([1.3476081611942619, 1.0569015062308722]),
array([1.193072451356013, 1.2032541581093719]),
array([1.4714995248039464, 0.8615634661755818, 1.1440298854625208, 1.490642709740996, 0.9611934454758262, 1.4289629227856673, 1.3607700319141083, 1.2882426953033919, 1.2107926092457209, 1.5001722756438354, 0.9113408338265481, 0.9375621256352007, 1.3850105138155229, 0.806832047442703, 1.5987768409925212, 1.5438051435409001]),
array([0.9328484449781447, 1.2725867441520795, 0.07386663142680619, 0.0]),
array([1.05719833470642, 1.1221840628228796, 1.0151743318755748, 1.1038125233716953, 0.3962172596348543, 0.0]),
array([0.9038059624846572, 1.117883182876222, 0.9719919038931799, 0.819799510189932, 0.6305822563688234, 0.5481427244599344, 0.0]),
array([0.8850508605122346, 0.9412194704681798, 0.6827149730875203, 0.6927836398083123, 0.0]),
array([0.40801233253294417, 0.4276528696878785, 0.0]),
array([0.6478887279527417, 0.5862076071704851, 0.0]),
array([0.12822166650352185, 0.0]),
array([0.238562813538581, 0.0])
]
d = [data_1]
names = ["7"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T3', 'T5', 'T6', 'T7', 'T8', 'T10', 'T11', 'T13', 'T14', 'T16', 'T17', 'T19', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T46', 'T47', 'T49', 'T50', 'T51', 'T52', 'T54', 'T56', 'T57', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T74', 'T75', 'T76', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T87', 'T88', 'T89', 'T92', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T102', 'T103', 'T104', 'T105', 'T106', 'T107', 'T108', 'T109', 'T110', 'T111', 'T113', 'T114', 'T115', 'T117', 'T119', 'T120', 'T121', 'T122', 'T123', 'T125', 'T127', 'T128', 'T129', 'T130', 'T131', 'T133', 'T134', 'T135', 'T136', 'T140', 'T141', 'T142', 'T143', 'T145', 'T146', 'T147', 'T149', 'T150', 'T151', 'T153', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T163', 'T164', 'T166', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T177', 'T178', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T200', 'T201', 'T202', 'T203', 'T204', 'T205', 'T209', 'T211']
def get_taxa_names(): return taxa_names