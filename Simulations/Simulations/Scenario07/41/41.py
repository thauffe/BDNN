#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.8466504213463]),
array([34.104381886090216, 33.2967054358159]),
array([33.33251390846006]),
array([30.16280644141729]),
array([31.1847463582032]),
array([30.631456405146412, 32.412225873735885, 30.444696637752262, 27.687920099289066]),
array([30.147162345777634, 30.267649075032647, 26.304812410589335, 22.727646233374383, 20.71483925131941, 22.808158960712753, 22.00591402513972, 22.07405988539186, 20.312768048529303, 16.542257821319595, 17.13692861748759, 12.573127125500934, 12.608204060592909, 7.682538103997176, 10.528828431283618, 6.624311494453509, 6.217797095547659, 6.877620462872464, 4.837435321099883, 4.917065346879349, 4.828942163480342, 4.811050985821149, 0.0]),
array([30.698096021568592, 29.579637180416583, 29.914321994726134, 29.97557638875234, 29.84636947493868]),
array([28.208499267478373]),
array([19.248936171314533, 0.0]),
array([28.41279091986776]),
array([29.311197641850473, 28.72188631907024]),
array([28.24705746285551, 24.254504766267363, 23.94456816090055]),
array([24.659542562050603, 24.64133440203319, 26.116650331319455, 20.821686870635535, 21.61894642474689, 21.709531996632013, 22.718954499902363, 22.593408043351346, 21.36149649775642, 19.792282265624948]),
array([24.35406666322498, 21.46078790055902, 21.198166378289827, 20.804113888999595]),
array([24.193176641765994, 22.888156649208703]),
array([25.388411789987224]),
array([23.434805725863445, 25.86465970098563, 22.140313140436145, 23.00428454698666, 21.69514616375579, 22.579463319185077]),
array([24.66921140358251, 24.396624618678985, 23.892439231022777, 23.344561281693288, 21.26125048409106, 21.786840658518166, 22.36067752940285, 21.206975020200694, 21.67774796893739, 21.103086612219233, 22.21581521653521, 22.23732353279272, 21.16809977375918, 21.670415566318212, 22.240995968854193, 22.628947147675067]),
array([21.529751926640472, 20.742001273250466]),
array([24.39400810560032]),
array([23.704904600971755]),
array([22.59353885894953, 21.27959738406603]),
array([22.89614496744572, 22.549908292147794, 22.05856920400879, 20.482889654319766, 21.88603602839586, 20.821302356924843, 20.928959102854627, 21.98420754485475, 22.65526196164658, 20.629607148338433, 22.249803782798416, 21.950944493139513, 22.081789894152745, 21.959964916119556, 20.813587721639188, 21.753728729833274, 20.79974571856053, 20.348353398387275]),
array([21.56708302279155]),
array([22.359743349725594, 21.974495145679352, 22.27135751097031]),
array([21.367250180778708, 20.612219783285305, 21.707082197821595, 21.481561896258206, 20.62447201221322, 18.067219554054056, 18.793364389920765, 13.939029668440703, 15.902833162485779, 14.437238663725505]),
array([21.533763434498546, 21.898200951146702]),
array([21.837266628383848]),
array([20.873758732480276, 20.6581788782872]),
array([20.46657196085554]),
array([20.550777054178564, 19.698397820496226, 12.881337280226553]),
array([20.601921420018286, 20.39608981405966, 16.686902784981584, 16.499666976418727, 15.152503213621102, 15.694378414616093, 12.494624965137957, 12.560163087631846, 12.17133917546149, 13.080759651833981, 12.241821554490793, 13.691134648803388, 13.036362864250298, 12.83290472463306, 12.885987581113556, 8.919395405179039, 7.286285580488245, 8.01587578915024, 7.501398742033251, 11.57942270308232, 10.380154764867045, 9.706134691151128, 10.55134854500838, 11.099629412308033, 6.2766426515405485, 5.851524361951184, 6.542582450976766, 6.09245207913329, 6.324148107938846, 6.549029037702607, 6.1871826395323275, 4.71595997815406, 5.051223355966001, 4.16954450055468, 4.411650632545975, 3.9188329652666125, 5.299953388603737, 3.9691536324325867, 5.095021249691301, 3.8307775964010533, 4.405124907549211, 4.688818768163138, 3.1612185155665284, 3.4709303023050158, 1.8459380757219552, 0.5844983731870597, 0.0]),
array([18.890123174757804]),
array([19.576515867647245, 18.208212581739247, 17.73705534978257, 15.477583039316277]),
array([19.160913000348817, 15.945919996986795, 13.443782530000176, 13.51671170895093, 13.504945880706591]),
array([17.857636403021644, 14.007246311288236, 11.945081143137783, 13.046045310657872, 10.769466903279918]),
array([13.934519113814595]),
array([16.80955873554663, 15.646067527580646, 13.654189279380102, 11.679309563236517, 9.108686531279737, 9.26971011909167, 5.972152197867211, 6.202041423575508, 3.787330276156795, 3.7201205452608335, 5.309721536944007, 5.031885128738435, 4.2226717953019595, 4.413089748783063]),
array([17.377363383135112]),
array([16.905791928836358]),
array([14.01820961019991, 13.214839792067044]),
array([15.90999821110283]),
array([15.429603672797496]),
array([12.545401253736397, 12.63149044485967, 13.537109226602114, 12.767344117637911, 13.416110655285253, 11.912409182045074, 12.500367333280279, 10.232132850576738, 11.05752858267581, 10.970153787450839, 9.53225677105994, 10.242073327836465]),
array([13.836361910627271, 14.05761844411692, 13.368719442649263]),
array([14.619572287363503, 12.66615036394896, 12.85444473234627, 13.119993394954063, 13.155502193654767]),
array([14.727028702829273, 14.766801710029887, 13.290428520271806, 13.063181821349147, 13.098236600786128]),
array([13.684561281004092, 12.24359047791968, 13.690093581276253, 12.03500505989322, 13.778582138051238, 9.636728540051712, 11.046134194717098, 9.200084958749338, 8.513268094514883, 8.986176421365176, 9.14476560517429, 8.087051394647842, 5.651605480900409, 6.271494061319875, 5.655537040957404, 5.4954187538387735, 7.042927086568869, 6.677464877443896, 6.735204463838901, 6.129098164612247, 6.779304491607713, 5.9962367802475445, 7.153526800799061, 6.237828224762306, 5.0630466114312656, 4.659953259694852, 4.920798169188296, 3.696834460403685, 3.8719629211276247, 4.045814891907933, 4.360175002710578, 4.520747859820105, 4.192828117314136]),
array([13.52581354215017, 12.963480213303484, 11.984958403576682, 13.491492613950173, 11.98247871064719, 12.232020749811078, 11.66769476855007, 12.19919131116873, 9.707650466211934, 9.64638199907138]),
array([13.96340695431248]),
array([13.765613185277285, 11.649890654792232, 13.488737292108544, 12.317765753885416, 12.768555990387167, 10.670353236590172, 10.980771209646457, 10.921975034383804, 11.452079912982164]),
array([11.880900925198844, 12.020691714313731, 9.493130828017422]),
array([13.224343634348697, 12.870398120184296, 11.954296448044804, 13.515059298387742]),
array([12.168044882407735, 11.819767770067202, 12.03101550576968, 12.960279203832977, 13.27302474842987, 12.091320429880291, 12.480843444338097, 12.94211360770604, 9.98370801482205, 11.420385748346522, 10.253548563642914, 11.429239703779855, 9.394161926234464, 10.727737332985159, 10.33216626760368]),
array([9.265759109487888, 5.96864273679509, 7.032287773370502, 5.31397654263576, 4.665714080530319, 4.232045745516049, 2.3848026763185723]),
array([3.6675009807809875, 4.781145688181448, 3.9658422319618536, 3.851847258332213, 1.6505122260242906]),
array([11.67391126494992, 10.833392822943573, 9.731304364496726, 9.499672622661482]),
array([12.091337508450714, 12.651458801461565, 12.65997314483602, 12.571237638029128]),
array([11.945526331781684, 10.73731369448111, 7.752189379881892, 7.671001307817347, 11.57423864928328, 10.512447642785197]),
array([11.854941073034825, 8.734761458504508, 9.406979185144587, 10.494712637332306, 5.425233198948656, 5.6442506243179364, 6.295831183831441, 4.049060662369471, 4.084999539329881, 3.8912378239750813, 3.9360385831967832, 4.852102992653005, 4.303713590365026, 4.961815317338183, 0.8405689029820449, 0.6123610515411168, 0.0943896780687436, 0.0]),
array([8.513258731965896]),
array([10.005002252129044]),
array([8.507372646618428, 11.12933240924693, 10.547379113369951, 7.102727295232575, 5.1013366146992345, 3.068827599643919]),
array([11.48638821937229, 8.740691407161787, 10.989381213703899, 8.779074114395563, 11.169637758193353, 9.564505166396964, 5.752843935950114, 5.758132436632136, 5.389905282445007, 6.010365323888519, 5.17068168795157, 4.6144474839026515, 3.615864542981334, 2.1782495161942563, 1.9504842590081992, 0.0]),
array([10.644260785395911]),
array([10.195296943886039, 9.315374368975736, 8.506119422986371, 11.082507605273587, 9.602886860965768, 7.9114879763659065, 10.046049595854841, 5.659952190828465, 6.621893450800158, 6.781721841745308, 6.238157045445428, 7.0440094401198, 4.708707902231179, 5.190857414156947, 4.870053678339824, 4.321563033944664, 4.052837217494908, 4.919612520341287, 3.1966510553077185, 3.301941666972327, 2.004373562742385]),
array([10.866511571846692, 9.365405845958863, 7.451065834832455, 10.364120932799093, 8.270028825931846, 10.40062972796511, 7.31258801970309, 9.149255914348762, 11.144147189849562, 7.580029387197097, 6.217040165958042, 6.363214234517085, 5.507506335309158, 5.518281775438592, 6.756478519575145, 6.691079677398642, 6.50892068202801, 5.438388351255626, 3.945739497025622, 5.127881926001613, 4.373772738998789, 3.949646764523934, 5.225507218414212, 3.967010756649346, 5.2742367684030595, 5.326347713893041, 4.694462748428787, 5.237929554226069, 4.030237668112546, 4.516896405791059, 4.216571781746687, 4.747875072509471, 4.2232076047938065]),
array([9.717741157897667]),
array([9.327783752442798, 6.588581795775764, 6.15708843484966]),
array([8.6725843922232, 10.321402547190573, 10.549108054115417, 9.19969418992311, 5.7559324192431145, 7.108244170410842, 4.364087643575813, 4.206738484524928, 4.468920200240255, 4.837468294447694, 4.277403575826413, 4.123638065170122, 5.097550997501475, 4.876854921632051, 5.0973801746800405, 4.764647628444856, 4.132838688500854, 4.1441793476959115, 3.133126492193533, 2.6534124662345695, 2.227761185702817, 1.9278951464276182, 0.0]),
array([8.122693362391072, 7.643302991643006, 10.351424787700221]),
array([8.526152221167324, 10.090758044599585, 8.97289585661147, 8.97080170258383]),
array([8.491547823015758, 7.546285271420066, 8.448236398106816, 7.191870925754278]),
array([10.370106527499997, 7.5162482507992365, 8.06863304104143]),
array([9.811034616879002, 8.767557756332998, 8.559331478218086, 8.175467323380364, 7.132577960719236, 6.4011130008471255, 6.3128922690846725, 5.772423924229071, 6.553561769472092, 6.122812477425077, 6.139147599260907, 6.257737262708309, 5.254646255795642, 4.543665109190582, 4.952583182357779, 4.8250797208674125, 4.418852918249401, 4.688832700059461, 5.109200613687115, 5.141488026563267, 4.727940523351918, 4.154999529793593, 4.942484408662991, 4.287942292864967, 5.235686257793854, 5.161230522359507]),
array([9.202710823660773, 8.708533980471268, 7.797541427989483, 6.358645940307706, 5.409398045851296, 5.78452437517949, 6.768410568404634, 5.490323877345913, 5.568641877964184, 5.338585011260229, 4.4519329284248546, 4.314569214125218, 4.708364972541715, 4.553210052606103, 4.277715247783025, 3.9460562865258058, 2.5966059784307642, 3.169401929951776, 2.3329519781786257, 2.0748484703356223]),
array([8.810024922984065, 7.750944939315374, 10.175488249029685, 8.000530838309125, 9.695737306951779, 6.552313659128864, 5.53356071611822, 6.162579102336998, 5.856659735820076, 6.382575732208605, 6.0909552074390945, 7.203765087347207, 3.6787851784154517, 4.750655047981429, 4.839090159252249, 4.108250970421399, 3.8576215645750356, 4.185410587700666, 4.946743792373203, 3.9374570378538096]),
array([10.090991255119723, 8.170120108488561, 9.798739986406854, 5.895057430580588, 5.44511745480828, 5.427900218713908, 5.513577102910813, 6.0732526293971745, 5.576068829406637, 7.216204612516063, 6.840567754040564, 4.800919605599109, 5.112152560634091, 5.041218524297688]),
array([9.671317827835182]),
array([6.9339737209689885, 6.051708662286292, 6.999497634055061, 6.889054955516533, 5.279907777074399, 4.666920375024004, 4.66253558962929]),
array([8.123932763777503, 9.138389305626262, 9.585331594104028, 7.564835305830723, 5.9779643036572825, 6.720227440152046, 7.053637737820919, 6.336938112891658, 6.774280399170422, 7.2290196852538, 5.7404583386642205, 5.280246100658904]),
array([9.176975499574427, 7.906448923703834, 8.324505540618906, 5.981557931213609, 5.358815549284338, 6.285626019416187, 5.812799650119162, 6.437116465946723, 5.5332265242990495, 6.636750684146478, 5.263057856035801, 5.313729906841062]),
array([6.429834771466491, 5.4411090018812125, 5.261529134630503, 4.528391633179294, 3.5461244705306516, 3.1370037779252495, 0.9494992666044061, 0.0]),
array([8.798413479734457, 6.328339054085019, 5.6310067600316875, 3.712591224578519, 4.907675808765929, 4.850793182628063, 4.472425718797553, 4.912149019839396]),
array([5.41603504715747, 4.8544326114974385, 4.6820388596584674, 4.75278004969817, 2.155811498427243, 0.0]),
array([6.193296446786372, 6.916740349741569, 6.834378766936074, 5.988845793442105, 4.832191164188789, 3.619259213364394, 4.533056726438049, 4.491658749723933, 3.1619733712612805, 1.959019961803027]),
array([8.779341644931224, 7.0060258406440035]),
array([8.808502045597805, 6.72291066321759]),
array([7.695748110750143, 7.888569036871419, 7.443263726411656, 7.41848805913535, 7.07967698843121, 6.235819474819705, 4.7392600393207065, 4.719925493426599, 4.79000411864148, 3.673572260200703, 3.946780486999092, 3.615992762412933, 4.271060103299572, 3.053304474295985, 3.4782852081253397]),
array([8.031456681733546, 6.143103547751429, 6.974092787940731, 7.164660129994905]),
array([6.4535410709463665, 5.442643520707748, 6.904677166957226, 5.637275261022127, 6.929789524507204, 5.632392839752038, 5.087382588634244, 3.7829766608163125, 4.687746961744523, 3.757995318052448, 4.160390773890637, 3.976064793352976, 4.916766870137065, 5.081478497625529, 3.9531520722693143, 4.973403821012274, 4.219827636173739, 3.12155156654603, 2.7765239060573, 2.1813041625895666, 1.9199421911965926, 2.2121873217294437, 2.3906982632485545, 1.9754861320228005]),
array([7.671210599580079, 7.469321372488425, 6.946946591978661]),
array([7.6190921284236275, 5.826455162928271, 0.0]),
array([6.621019765201184]),
array([5.356444092683351, 5.701533693131373, 7.117408322048137, 5.056315974161081, 5.241766192920433, 5.252548968545943]),
array([3.85373858403573, 3.6071597031815568, 3.773353729185255, 3.924037259805914, 3.8119823397997807, 3.440454370469027, 1.8655679612298393]),
array([5.940768324547471, 4.008683133347446, 4.476021749572146, 3.44391409884988, 2.37036618957339, 0.6025300485836474, 0.0]),
array([5.586244886109488, 5.50168821558556, 6.027948411876867, 4.093730954747306, 4.380178958943331, 4.51494491275756, 4.05787142063151, 5.296908585136363, 5.317894985930343, 4.811170985242536, 2.113697016695271, 1.9052633538440993, 2.084288686198187, 1.5737032578030936, 0.2049495719418234, 0.0517760150849927, 0.0]),
array([6.4785569814632336, 4.903271103700335, 4.694545005088308, 4.218411037265788, 4.071509019686581, 5.0022937928218605, 3.190992208985172, 3.079203159295394, 2.237512551511147, 0.648427408535277, 0.0]),
array([5.587154886663953, 6.410884325885436, 5.106221740565763, 4.309878610263455, 4.229459073943758, 3.749413689265661, 3.057128179262608]),
array([5.499812853020411, 4.516900480595175, 3.7438263571342856]),
array([5.926596224337456, 6.60257060482093, 5.95729294284274, 6.489024603947899, 5.622411567434934, 5.934979321926184, 6.293185705892956, 5.823474376282155, 6.02430583808654, 6.467567214638631, 6.41183017951374, 6.513433341586618, 5.295053135519042, 4.558270743637189, 4.995421675325324, 4.43072353582502, 4.260145298349697, 4.493999345789389, 4.263368193089289, 4.2356477303713485, 4.806639706302604, 4.769871055759654, 4.239622762572104, 4.864940714201316, 4.343214566484797]),
array([6.190128944802476, 6.63596132237444]),
array([5.577444288734809, 6.464536894599792, 5.695417529388803, 5.516336357956435, 6.638250670821783, 5.413494080916042, 5.575783445329195, 5.1275305894397345, 5.1752784462560735, 3.606843040923123, 4.830380246385137, 4.528623140408498, 4.939932653501101, 4.465259171870575, 5.225315026618854, 4.575702944594168, 2.9099067535593326, 1.0461672173836618, 0.0]),
array([5.646044340854132, 5.724518790487737, 4.458121413720406, 4.88753643027621, 4.4866322427160235, 3.8629532206200032, 4.407172468840754, 3.6919687193837247, 4.9352546198054945, 3.845452035300939, 3.2782627669098883]),
array([5.1313081204537, 5.221784286962428, 4.902332529229714, 0.011396065262457386, 0.0]),
array([5.978860553135463, 5.34962935531411, 4.513420325987161, 4.28473454980449, 4.205296066033252, 4.914451014352578]),
array([5.918501817095272, 6.134610466780179, 6.206454532109891, 4.841610588099327, 4.088926942893046, 5.0432375834707734, 4.336013789247177, 4.475913500492945, 4.0615480174038066, 4.520239682715586, 3.836013400698982, 3.745181766447905, 4.3000930289082335, 4.0593959056481745, 3.6203100764566956, 4.541467216693903, 4.125360986023825, 3.1490884542672344, 2.835563311197901, 2.771811705521137, 2.164636057608337, 2.1603201785677935, 1.0241494517271228, 1.4382254802061791, 1.2552324721078434]),
array([5.728911586937061, 5.462599259446887, 5.7506970381015465, 4.465345595397226, 3.8266271396300144, 3.9820722595993225, 3.661984537909671, 3.9241041391009768, 4.6653779618643085, 3.662268934346665, 4.786118401885002, 4.1813115957455445]),
array([5.760208362919934, 5.463195618312713, 5.434904632626831, 5.543735916936865, 5.3731511492123385, 5.964900440060189, 5.198921262903702, 5.290687587429827, 4.61931148588643, 5.274307807647933, 4.241382256833113, 4.1423149815176545, 4.729104849830833]),
array([4.082830029980679, 4.0938986013141445, 4.727458739365639]),
array([5.392061055889114, 5.318443962055674]),
array([5.613024877242612, 5.494880984771376, 3.724690662528743, 5.128759098018175, 3.7921502002209673, 4.877573407231366, 4.728853513350088, 3.8271537330801504, 4.530081619199571, 2.6316715178314745, 2.2679897792835657]),
array([5.508018983738559, 5.465329847155701, 4.209398334765954, 5.223144240212699, 4.397917637644681, 3.9953206237135523, 5.120625576414456, 3.6214825702736473, 4.21269084620663, 4.129115632666485, 4.679811942790794, 4.760826334281098, 4.9624968948991315, 4.345209053466545, 3.3508864385209347, 3.1422416043070815, 2.4173948044517686]),
array([5.7114292575806385]),
array([4.2726208524355025, 3.774271325451235, 5.252638612686182, 3.8027553481410825, 4.295672498365988]),
array([4.968406726553047, 5.244529330115692, 4.216235356951112, 0.0]),
array([4.7529977558018235, 4.076652393519934, 4.708500854275201, 1.9140162936807645, 0.436262497162558, 0.4716746547367983, 0.0]),
array([4.830415497366034, 4.914579805517404]),
array([4.024818530412787, 3.908678398130175, 3.74456581566968, 4.300518245173458]),
array([4.503917268344954, 3.910248528861919, 4.184519380975189, 4.713198106168055]),
array([3.7990693291619504, 4.270830752903565, 4.78785690511765, 4.689858093361038, 4.7012983517124995, 3.9256660804928036]),
array([4.77333302657124, 3.8187316081533007, 3.38079211277285, 3.3662438936041275, 3.417534180580615, 2.493043799031946, 2.2710654356037687, 1.2747335630883714]),
array([4.675060775009926, 4.702180355097814, 4.722854542636836]),
array([3.9731100448926036, 4.08625973589825, 4.65168856260286, 4.477166425889048, 4.077750286858974, 4.145328847944311, 4.49202363058336, 2.8194680285027602, 3.3844848462743724, 3.056502870105537, 2.222142945985147, 1.179297976415429, 0.0]),
array([4.484082431129321, 4.509897812508977, 4.464326990351622, 3.672114450741148, 3.915924154181591, 3.9452528905622897, 4.390110058055502, 3.767981654462658, 4.489050479847636, 3.8326109590051782, 2.6037155649530814, 3.504464321859582, 2.787102005117457, 3.4813334917810526, 3.364590425557358, 2.033058546785338, 2.318862119134307, 1.4348175045967437, 1.495399775699087, 0.9761318984521241, 1.2545102510077037, 1.61890770612366, 0.936290895743368, 0.0]),
array([4.028663657008165, 4.616204683062397, 3.619050496351746, 3.822500870959543, 2.455175619602361, 1.2817091175507191]),
array([3.6638332151579918, 3.953212077153771, 3.9779985085910097, 3.825784907408518, 3.626734947872281, 3.950159927966261, 3.7998646371077562, 3.886369164911434, 3.7534998128197405, 4.093920064364623, 3.783153942887169, 3.925167381077774, 3.7070390868644143, 4.011738227677237, 2.944342826804416, 2.9413565022112276, 3.0751071566562183, 3.5269381506638786, 3.182630465201225, 3.4751807227223854, 2.8457671993279443, 3.355202621949261]),
array([3.9873736572815255, 3.886611484183722, 4.116460397243626, 3.7238027132027995]),
array([2.6296681870503535, 2.9739820284689906, 2.69620384163917, 2.929649819757823, 3.4869045787869517, 2.3436302673165996, 2.35884146753329]),
array([3.92697863803539, 3.87775924641705, 3.9442571823089745, 3.34067983918516, 1.8334410562884627, 1.6338025521873727, 0.7325753161836196, 0.0]),
array([2.7300587654393103, 0.0]),
array([3.808492964351382, 3.7922293281337236, 3.5371666565727558, 3.1919672350740624, 2.1221439201205308, 2.1686074445400654, 2.206166325996851]),
array([1.8788259969968413, 0.0]),
array([2.894117287896854, 3.516625791904576, 2.2168462519997876, 1.931408555347343, 2.2282210856602074, 1.5933584698134924, 0.8117616455833144, 0.635707414479032, 0.0]),
array([3.2204222265138953, 3.1804814219101334, 2.9599364600905806, 2.0022935865413434, 2.428442725609201, 1.4900060817845555, 0.48552905497926807, 0.0]),
array([2.0375565752637566, 0.0]),
array([3.0095197573947163]),
array([2.842231350386936, 2.805086269019227, 2.9005674364812495, 2.998911283769565, 2.62947456775856, 2.541758013313505]),
array([2.7898862889944036, 2.989403616481339, 2.59981348497461, 2.991962470068607, 1.8941168445665606, 1.3892943356627256, 0.6265458210092133, 0.0]),
array([2.835320691977794, 3.037906460038047, 2.231161566282158, 1.8772457413764174, 0.579275635509796, 0.0]),
array([1.6757480093132133]),
array([0.16222963557362002, 0.0]),
array([2.8957913284877077]),
array([2.84206471793341, 2.0348664753320618, 1.9445248276000169]),
array([2.6752777507267016, 2.6783628591409174, 2.068139595940951]),
array([1.8643613901274323, 1.6587727400195627, 1.770712057547879, 0.9978164102073431, 0.19792684751810552, 0.7038089075530375, 0.31584234451671583, 0.15896874782778148, 0.0]),
array([2.5961422474076876, 1.2060469498243087, 1.527995847637774, 1.3444306280671467, 1.766494858940762, 0.4405765389618627]),
array([1.2386920322562494, 0.8471312967254756, 0.9050680449848152, 0.6018369273424452]),
array([1.4709241274757225, 0.59561069639399]),
array([1.9289068334465393, 1.6314313735748134, 1.7251579359954652]),
array([1.7866840877548351]),
array([0.9053228525001811, 0.0]),
array([1.4526108880300646, 1.40994166642606, 1.2922182885997846, 0.0]),
array([0.9817645795643485, 1.1710999623054832]),
array([1.6026472518387422, 1.7319115214574532, 1.244879916056868, 0.0]),
array([1.1025256673291142, 1.4160870746547898, 1.6155204421734, 1.2746820244154138, 1.5157845197054531, 0.14762792027566862, 0.31178577871121643, 0.0]),
array([1.5104717410482484, 0.5750999192973445, 0.0]),
array([0.23322808058268651]),
array([0.7703484352300258, 0.33872380756258125, 0.0]),
array([0.4859365123772158, 0.0]),
array([0.9537260924497591, 0.6038122572073736, 0.0]),
array([0.5424913281600434, 0.4147549930068938, 0.0]),
array([0.521014134639973, 0.036783862795905536, 0.0]),
array([0.2824285383654878, 0.0]),
array([0.6498462781183352, 0.0]),
array([0.28232769121788887, 0.0]),
array([0.15837478658702742, 0.0])
]
d = [data_1]
names = ["41"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T7', 'T8', 'T9', 'T12', 'T13', 'T15', 'T18', 'T19', 'T20', 'T22', 'T23', 'T26', 'T27', 'T28', 'T30', 'T31', 'T33', 'T35', 'T37', 'T38', 'T41', 'T42', 'T43', 'T45', 'T46', 'T47', 'T48', 'T49', 'T50', 'T52', 'T53', 'T55', 'T56', 'T58', 'T59', 'T60', 'T63', 'T64', 'T66', 'T67', 'T69', 'T70', 'T71', 'T72', 'T73', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T84', 'T86', 'T87', 'T88', 'T89', 'T90', 'T91', 'T92', 'T94', 'T95', 'T96', 'T99', 'T100', 'T101', 'T102', 'T103', 'T105', 'T106', 'T108', 'T109', 'T110', 'T111', 'T113', 'T114', 'T117', 'T118', 'T119', 'T121', 'T122', 'T123', 'T124', 'T127', 'T128', 'T129', 'T130', 'T131', 'T133', 'T136', 'T137', 'T138', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T179', 'T180', 'T182', 'T183', 'T185', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T197', 'T199', 'T200', 'T201', 'T203', 'T204', 'T206', 'T208', 'T211', 'T212', 'T213', 'T216', 'T219', 'T221', 'T222', 'T223', 'T225', 'T226', 'T228', 'T229', 'T233']
def get_taxa_names(): return taxa_names