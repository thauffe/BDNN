#!/usr/bin/env python
from numpy import *
data_1 = [
array([32.55542663219098, 33.03599989974175, 31.67799553238925, 28.14875808464607, 32.322418199422195, 32.174221174512404, 34.72283040998799, 33.96621646497143, 28.428230451553873, 30.55739983527973, 28.933872240965727, 31.292721734737896, 33.835204950721966, 29.61996236087057, 28.262719008630995, 31.23236369573901, 34.99906108045448, 30.19372083497166, 29.24478491012465, 32.88809508905408, 28.899564709601435, 33.796523439690766, 28.22153070381135, 33.42507682842427, 32.42396338786665, 33.80934845284733, 34.89078881335932, 32.23233969909524, 29.43119726708055, 29.59471752453777, 29.761714045100305, 32.5934171807726, 33.75224559471837, 34.65117691382243, 32.944227422956835, 33.22336547598867, 28.94291813808091, 29.142736753006396, 33.6929702663341, 31.674021450240293, 31.36533791828909, 29.145057246964733, 29.37248854611454, 33.142387007125855, 28.86500160897685, 31.647246962558135, 33.821108728084184, 34.94424912540354, 32.340598056389176, 28.57040833763894, 34.233273972162024, 33.38365931531666, 32.23615442906927, 29.026464810949648, 33.0097990125459, 34.80462415958771, 30.05339535774764, 28.397587595643422, 34.58201794971967, 29.325370425467867, 29.580577860209146, 31.13616261208433, 28.728932052634832, 30.911040369860622, 28.623175831042374, 30.51339684248801, 26.40103897810852, 26.40070115055565, 24.140349558612773, 23.925543502025395, 26.004489824015494, 25.203988582002996, 25.499402272949126, 22.31710581121743, 22.85224024677864, 20.499639640041348, 22.37955434559967, 22.07280590977023, 19.656685792581985, 18.686380944596262]),
array([28.63341660885594, 30.056128155064332, 29.60717941211724, 30.295345775075447, 28.24250548263479, 28.424936610588837, 28.180668992163604, 29.572763613792926, 29.138746375488097, 29.95209392271556, 30.719419268217095, 25.786076136046518, 23.403497131672996, 25.412364400311063, 27.004738549682816, 27.788932745287287, 22.146648028054408, 22.791358887107837, 20.666924828416395, 20.726101443428583, 22.04439304032115, 22.53897176893228, 16.335527329063623, 18.870464335393134, 19.591161607051404, 16.79832995974635, 13.200240946080646, 13.21702075820962, 13.378611328651125, 7.600312518379634, 9.45322752637325, 7.653536106748973, 8.168136771381167, 8.863735145301293, 8.087418739558169, 6.905058634629244, 5.7752242143227495, 5.655645712196375, 6.9447370465938265, 5.233126390188329, 3.7299897253447902, 5.248529373724975, 4.325222147296848, 4.867249770874851, 5.273255134318631, 4.4786406197553745, 4.143633536303877, 3.601836263972693, 3.572338043383816, 0.39651205973079306, 0.5631468498337808, 0.40596990322112547, 0.42461308189738656, 0.04678011057424329, 0.0]),
array([28.215935653116876, 28.801828395457285, 29.085720429973893, 28.267533252396355, 29.465692864924744, 28.446028928200857, 28.769319301401254, 27.763668279381683, 27.931737173258355]),
array([28.239841209209665, 29.190615301889782, 28.948718477210225, 28.703980760832998, 28.219132530795772, 28.845178240099504, 28.413037722836577, 28.661500432990948, 29.210851608752655, 28.634831341318158, 24.61396072656592, 25.116941857898873, 25.398588994526875, 26.313411051933862, 25.26817175595337, 24.185597204727813, 21.42403891247058, 21.087716420316205, 21.102536779452077, 22.294745551734298, 22.008312961185673, 20.952571195082, 22.23035318393315, 22.690208423174596, 20.30598442902202, 18.46341562691607, 20.147971798812602]),
array([28.73053277916882, 28.522184916970833, 28.705177601443406, 29.068643896473013, 28.445156076720853, 28.30973707435405, 28.717471878681607, 29.063863540834987, 29.125512519322132, 28.717771771858988, 28.883042631618345, 28.350321108601516, 28.441684988939443, 29.049105662177364, 29.213673882890056, 29.09851581909821, 27.387432590066396, 27.987677887520675]),
array([28.261963248838338, 28.523989532032775, 28.52699930640631, 28.91216952988962, 28.829033949153565, 28.590205005806148, 28.503967899387472, 28.3539012374805, 28.52617573452115, 28.59845310949401, 28.847751458220213, 28.247969159727035, 28.834537109846927, 26.350359364248867, 25.945078165095875, 27.626218250107147, 23.630524329344485, 23.9733212613235, 24.224804497698535, 21.20675214184396, 20.61855430243443, 22.564248424904218, 21.529074264208674, 20.48070974244156, 21.058303460808947, 22.189541713341637, 22.05590509839113, 20.652492685026292]),
array([27.56669514727031, 21.819513531416344, 21.27512255566834, 22.145132736815626, 22.575329452730365, 20.312246678090705, 16.95986811152543, 17.963783905596923, 13.535397097094284, 11.666702540196207, 12.129679023496973, 13.605872605712367, 12.763914057147012, 12.148825020732529, 13.465706415393797, 13.338925683138854, 13.801970928002628, 12.019655734384067, 13.19556797312747, 13.076186282160101, 12.866033228780445]),
array([25.142058766893168, 25.469110804343547, 25.706939919047247, 20.788722596913594, 20.797639939909683, 20.89054134716644, 21.33441934566417, 22.650588299815134, 21.12086761578935, 23.026124511208042, 21.679562293020997, 18.303153221967037, 17.03149628331611, 19.413771081864493, 17.140684357484616, 19.421523039267264, 15.49555669237947, 12.052913720970668, 11.657305618979196, 12.90603439143626, 12.553895027479388, 13.610655439108019, 12.81089780477048, 12.376709151313744, 11.740321445714988, 12.140535755537178, 12.368801347758389, 9.120431209345318, 9.67015124115949, 10.46213427508148, 11.340199053104845, 11.27408558089655, 5.679005957011093, 5.010712792617951, 5.248783764981539, 4.944305182727638, 5.078285805799166]),
array([26.613669429320268, 27.195584194778682, 27.506586499389073, 26.59799332543963, 27.26644372372555]),
array([27.213653717977113, 25.90851425967658, 27.200840403542273, 25.7028644005729, 22.514781095671307, 22.195731220912656, 21.89763720786022, 22.44230498073894, 21.779173990813597, 21.96655028005068, 20.764622467984463, 18.643103548287126]),
array([24.194051088765598]),
array([24.399886139653194, 26.098641050908466, 21.75412600241359, 21.451248424689787, 21.14855438290339, 21.289856174444054, 18.02618632761005, 20.072656133236183]),
array([23.833453580210197, 21.41221446837104, 22.45430129163078, 20.451214579628505, 19.71141213706171, 19.58511802520764]),
array([23.06054957769839, 22.093193071920584, 21.184982550200242, 22.353419942609207, 20.848576853177093, 21.826752280070266, 20.861127330183518, 21.046144568703813, 21.061082459936337, 21.846213932968713, 20.785342069859443, 22.574996222604142, 18.231988702251172, 17.22078738228731, 18.605192258618523, 19.145358640322048]),
array([23.019760394192957]),
array([21.817071563673423, 22.52038196284991, 21.791215472564122, 20.727917977218745, 22.257628503671647, 20.93311197333857, 22.73830292892627, 20.979413494324675, 20.946889919636227, 21.552020479899898, 21.514454042349406, 22.408644650829782, 21.917579975955285, 22.413248247056487, 20.471399167533317, 21.778470954497788, 19.24092715575833, 17.88237866629557, 16.81857752115814, 18.336105615231723, 16.487011259300683]),
array([23.39971810422452, 23.084998637328116, 20.94224768582056, 22.44143437379422, 20.787744642522114, 22.52477162841954, 19.632829395720364, 17.420583292144755, 14.736964097209256, 15.895611887691866, 14.207049282643792, 13.485409207989397]),
array([22.2583884864261, 22.148821972847045, 18.33689029136078, 18.856940608408774]),
array([21.34415364111757, 22.336632982650865, 23.010416179987327, 21.203382343287668, 17.59718133214544, 19.220144783903574, 20.132670684148433, 13.85762734850122, 11.810522304043793, 13.494210483993164, 12.06527746087355, 12.208445370107457, 11.257089458733772, 9.985000326670109, 9.878543041758967, 6.955158433904885, 4.404297943068829, 4.551950961934076, 4.752625447315032, 4.951764954527921]),
array([21.461636609383692, 22.074694050353038, 21.231667109823857, 21.54117065763271, 21.10076977237051, 22.020006765443092, 18.696742509303967, 19.064259739981434, 19.328172667160917, 19.484254772857593]),
array([22.19269203640216, 21.52095087877738, 22.33038165540677, 21.637624888820408, 21.406323283710314, 17.535995468468442, 14.47782714221912, 14.998228583760639, 14.791901490392227, 14.714073736742417]),
array([20.979678941234074]),
array([20.895518437927844, 20.477832528085294, 21.06259436169974, 20.939943408585332, 20.00800836070872, 13.377217380750032, 12.1550631750776, 12.517909886671989, 12.649964658063576, 12.861772544489725, 12.211836590428156, 12.219736529345703, 13.664114189315688, 13.227005350875137, 11.800364287658793, 13.336239819673128, 13.801558618347816, 10.292270553083796, 9.298184751221125, 11.48610168296891, 9.653226564859573, 5.931293597769004, 6.796247112379601, 4.61031863165529, 3.9884078154483156, 5.021060225500875, 4.31238759013902]),
array([20.89691913979615, 18.166072072225745, 20.4293138690846, 16.914531242345962, 16.371143630070634, 17.38606345182481, 18.523293106495583, 15.925549696845392, 11.919235611653088, 11.173120225710628, 7.606361457949992, 5.534345176314282, 6.778222373294586, 4.646061961086179, 4.422179809261124, 4.004677682412631, 3.3642756890837977, 2.9494152871164276, 2.47987887440338, 0.36766391612483207, 0.579647847597807, 0.44084034470363653, 0.0]),
array([20.57634522437643, 20.931359888794734, 21.009413386403516, 20.601330059360485, 20.050737275573198, 19.716261457799618]),
array([17.72496994871633]),
array([20.94492412622963, 12.931111827945038, 12.720338404420406, 12.056094064407597, 13.21652350470472, 13.427183188808721, 11.787996382313944, 12.338374157620574, 13.71899830628409, 11.889338223472082, 12.63298253459472, 9.030695808164701, 10.838571678921701]),
array([20.89196397055826, 20.487056639774003, 17.63687093768835, 17.464527801306414, 18.59460013229352, 19.946632616909444, 18.313172865693286, 16.38842455098689, 19.17726173715065, 16.394345275625703, 16.159145517523513, 18.039449006958677, 14.426571911944949, 14.328705666577184, 15.661677122179855, 14.103796893253604, 13.634086672322015, 12.972892140488593]),
array([17.925607469175635, 16.38805170872824, 17.741472189744364, 14.5532731598737, 14.595204341229644, 13.065048264475152, 13.492042992069365, 12.729456556737608, 13.29330426214098, 13.406898687094532, 11.689762968997115, 13.359898594161907, 13.642934301409856, 12.881181532840849, 11.922470851195046, 10.440251957569549, 4.399441197906401, 3.7453561888263094, 4.925382364276089, 5.048108266596048, 4.382806125937681, 0.4909045856666877, 0.6120220232818268, 0.0]),
array([13.192197043392426]),
array([17.153523732358302, 18.447868394589776, 12.237049445019654, 4.945333412194047, 0.4099824505791917, 0.0]),
array([16.522434327154983, 16.26113207518956, 16.686333525617997, 15.026530710014098, 15.316901753384814]),
array([16.05755573755212, 13.47053668872837, 7.5116571663800205, 4.454569055436208, 5.31981282878724, 3.8107580304938296, 4.492195399298923, 4.139769919751551, 1.903583658439887, 0.0]),
array([17.081957309617195, 14.320557137171136, 13.96284821584871, 14.818940150331649, 12.461190727663215, 12.4502564606777, 12.527987245385884, 11.700659164772183, 10.122382507408899, 10.157272992286094, 10.137740391137527, 10.303053001447934, 11.355512225172761, 9.696466834804843, 8.581304308921155, 9.81536350215389]),
array([13.920806134799784, 15.130744309983529, 15.256735333597328, 13.180849221358768, 12.702650127266864, 13.341477281942367, 12.631635918411067, 7.903461789223375, 10.474099033179337, 5.887097497640948, 6.937248986690982]),
array([16.177779134494838, 16.517778311763255, 14.56792416263822, 14.370555373141535, 13.92015207012066, 15.436497364232018, 13.127423570810533, 12.092504735558865, 12.453708000435727, 12.245612785331499, 13.509302047304653, 13.568902505980919, 11.726973785972698, 13.102482225721316, 12.339698194835751, 11.677530453452666, 13.723315978887166, 11.940635674220482, 13.005496382935798, 8.240105124001246, 9.450277791456365, 8.129267607138747, 10.69573748616752, 10.339095403190619, 9.946762435151271, 7.89454978559767, 7.4750353829299545, 9.618925959690408, 10.63799388048363, 7.478387487529094, 11.146826109230066, 6.179372255187941, 6.703144597113184, 6.779441741456139, 5.987847179565907, 5.744212547555002, 6.545321495843303, 5.960087703545467, 5.149251452395342, 5.021529918977679, 3.8363435349083392, 5.069452364811563, 4.005329364677339, 4.838428360913957, 4.675968529731152, 4.226135913624644, 3.8471257263304484, 4.738244496190898, 4.975166401945657, 3.8356954455198116, 3.9213566517796608, 3.806790508634167, 5.189314881677508, 2.7883902586364293, 2.4729595693643462, 2.296044009081365, 0.6898624369306813, 0.22088177210380033, 0.6114002577837003, 0.49874823574422095, 0.2533518384496031, 0.5442196890392055, 0.0]),
array([15.682182053574403, 14.520775817157476, 14.512337003417501, 15.184161829345427, 15.144502386317328, 14.678503398563421, 14.992533285013687, 14.954534972641246, 15.800177922275125, 14.225223673684487, 14.536989934472464, 13.678318253908813, 12.386096335155313, 12.410004673708526, 12.31884577830888, 13.137109712446524, 12.394072088157188, 12.840478059380574, 13.451801164541788, 13.79969268804605, 13.350721671249524, 13.011845619725003, 13.238916490033512]),
array([13.034182050249292, 13.06597907153526, 12.521035520176614, 13.331030434893192]),
array([16.22822253535485, 16.186223828480703, 15.134166039079734, 13.989330375107917, 13.848998479811073, 14.96852396289448, 15.057609153713654, 15.235636662509965, 15.201313527814863, 13.763672555755635, 13.770758555860507]),
array([15.987723114851683]),
array([11.804528151431182, 12.127282945384142, 11.957533832065678, 12.862419309976172, 12.971628283585789, 12.472182199691598, 13.76876582662299, 13.406866586345755, 11.685385863891536, 12.470394851201627, 12.946169080528424, 9.156905093218876, 7.21144617115475, 5.962035743701337, 6.986399328285952, 6.221315664067428, 5.745046455547996, 7.233871525903147, 4.111185937465208, 4.495649788746412, 4.668854820803655, 4.254531930408191, 4.942876653186232, 5.202759095654908, 2.656479452232256, 1.3440595684115415, 0.6764525682978902, 0.45460235015431383, 0.662158509086209, 0.47280798102431, 0.0]),
array([14.505881308989174, 14.529685955775255, 15.602044529421669, 13.631831211117298, 12.544092862103806, 13.190333308394468, 12.634680583772438, 13.684958153368122, 12.496197659074156, 11.967075581927839, 11.378012098704586]),
array([14.856885235391587, 15.706382420172288, 15.106991755946845, 13.845826043716478, 14.388927877148866, 13.768597965263886, 13.420480592614908, 13.482011744415688, 13.61769132355107, 13.561370738168417]),
array([12.388134204023313, 12.011424734685146]),
array([13.854872450588772, 13.890612813129636, 15.030257913524705, 13.029410563954942, 12.59536610400318, 7.566015868322747, 11.495870148257032, 6.2409618010565024, 5.873087181914027, 6.238887844632778, 5.884813572413297, 4.8599957847242, 4.140906642463275, 5.223176514431062, 1.0428173205973779, 0.6052044790507869, 0.0]),
array([14.516967886679744, 12.185430674073498, 12.300445055145955, 13.108986319076982, 13.06564640616878, 12.052258617870871, 13.56030672414889, 13.54946250293074, 13.00206782959485, 13.819104354306129, 12.996333013255354, 13.099342063527184, 12.09972369480825, 11.406002834012288, 11.15492109580609, 11.522455931879914]),
array([13.05070445380295]),
array([12.371966168716762, 13.785772036946819, 13.631284786079762, 11.842422100998167, 13.753060064355997, 13.561822925451825, 12.431813943743286, 13.596874214170137, 12.969362188317398, 13.315313576960383, 13.017380318197652, 12.16666603547637, 12.99963477215902, 13.772052157087257, 11.748944223572636, 12.075617605778397, 13.069068838193179, 12.74230299552643, 11.839772364952491, 13.733956740988019, 7.342239652878307, 11.40612641344124, 10.617175621664162, 10.799062275442072, 11.263615723583847, 6.819471020381329, 6.70949549744763, 6.0064696235361295, 7.155221734455573, 6.207109890202002, 4.421317756234752, 4.574314761017341, 3.8191824074218657, 3.611962999363014, 4.998961416974125, 5.134092129218211, 4.307829535510791, 3.9111592745212462, 2.7102738429216693, 2.4647085603394645, 2.1505060827859275, 1.1196852394714791, 0.26509866263016946, 0.7297358600393036, 0.1995465840498275, 0.08817174708797315, 0.0]),
array([11.716287943028622, 12.385706187277318, 10.472463713278982, 9.367026462804247, 7.160769309764514, 4.92771315001205, 3.1708010467166607, 2.820750447605178, 0.3261188446905401, 0.7649599651614181, 0.0]),
array([12.184579852940804, 12.211925895154828, 13.04775198160409, 13.38201957550298, 11.642850161946054, 7.601605845159467, 11.59798152503999, 7.716758476087144, 9.69547873852138, 7.944385257964421, 6.79618370610856, 6.219743620791381, 6.2663531464070745, 6.75552734415382, 5.461175799905146, 6.16833688419195, 5.0788549205242575, 5.233366020773989, 5.171501644030849]),
array([11.756976712644896, 13.568241752489882, 13.116919477093886, 13.556980762972382, 13.011080805154826, 10.267856243276386, 7.630736594968182, 5.646542651627597, 5.538424395964505, 6.216244819694592, 4.033428000653673, 4.8041959830499135, 4.761630621080144, 4.522240768059426, 4.101128841028483, 2.9351773643472048, 3.4053483547706542, 2.5773904341705203, 0.6238437407947073, 0.0]),
array([12.942597583368254, 12.322274472960578, 12.842594935887236]),
array([12.938879266336883, 12.255908620031832, 11.97132256246578, 13.038797918407175, 13.023037562873572, 12.647013574869634, 8.780966332941686, 10.679361956309908, 11.090094386292487, 11.237638690151991]),
array([12.089045593862396, 11.797871788704247, 12.698020128469409, 11.679365535896707, 12.815278895334814, 12.139350780116809, 11.68936234776213, 11.34765628827697, 9.371268823785368, 5.901465848139197, 5.6192623486397535, 6.325303532979567, 5.171470097119474, 5.27405863617491, 4.600845451325958, 3.390578761267387, 1.9216005261410953, 0.0]),
array([12.040436259728583, 11.648916733891003]),
array([12.210861740530747, 12.297447321215236, 8.628353748150367, 9.530730442790952, 9.804547833195306, 10.350625284819206, 7.995413137863434]),
array([11.705107220475483, 11.147262965934715, 9.403219201811376, 8.5924452053661]),
array([11.655161455311687, 11.634249956488796, 11.211629772610305, 7.00207432088177]),
array([11.698739912859727, 11.699187086392318, 11.954536834451742, 12.12393901331025, 12.172546057700757, 7.433445148774767, 9.65683163105594, 11.18506763114876, 9.243540767711366, 10.654497398216611, 7.8096702856699345, 7.1638920960389765, 5.4283916735952475, 5.3707054102109595, 4.305104735132426, 4.376471979440359, 4.401183071412468, 4.1761730572561575, 4.45981470063401, 4.513259411109978, 4.231628543740162, 4.252746078299953, 4.170829415847263, 4.93853076135476, 4.426695605679971, 5.205844392738364, 3.6644396003124156, 3.4099275058236467, 3.262838689780562, 1.9511818761560913, 2.3426019623931005, 1.7973823701298828, 1.5150482056882952, 0.878748209605978, 1.0557215654528334, 0.5397820279705827, 0.4867041070880458, 0.0]),
array([11.83680211162623, 12.282101969321815, 11.984102580535247, 10.803914698978101, 8.74031180125047, 6.949394487705171, 6.512257332106085, 4.480874142401324, 5.222312333902803]),
array([11.927876462951202, 12.179576874822699, 9.043242929683238, 9.549864903417962, 10.7086957427263, 10.763959600314692, 9.266904275089631, 3.80818820154902, 3.7443438673905414, 3.8240405702793066, 4.443191123427581, 4.349098480155537, 3.262287211432948, 0.5874539524307484]),
array([12.00728200568499, 11.784167729030791, 12.087452779468542, 10.401545835673284, 9.208450652895582]),
array([11.149391434354373, 8.595054173786465]),
array([11.646467741981866]),
array([7.9997355901466545, 8.53071843309953, 10.047259982475367, 9.842819170731252, 10.596133172675161, 10.443303890083145, 10.923822997749916, 3.810378941206716, 5.183116355275393, 1.3364335881949239, 0.21463481171796361, 0.4449132288786511, 0.0]),
array([11.255755022567156]),
array([10.985941968166992, 10.848965318474642]),
array([6.0629397810697725, 3.774542806885863, 5.03587239517524]),
array([7.549857838706504, 7.960023992294856, 6.536651316840869, 5.312223625332696, 5.197948378325103]),
array([11.075244738378558, 6.587033705975243, 5.01860159616527, 4.560976412205641, 4.884897581042986, 4.620429893534499, 4.395664066347481]),
array([10.884946761534938, 6.329083960363384, 4.923058492106815, 4.391350933716278, 1.5404397584933247]),
array([7.410461936383754, 7.47678895535349, 4.8388093849427305, 3.60802908507175, 3.9796695152317447, 4.632569255505938]),
array([7.358221320456698, 8.839281712563473, 10.174155079613799, 9.236765701098347, 10.241406665518907, 9.343098949199453, 10.491602099604721, 9.083838721015276]),
array([9.948854460956099, 10.278543449468257, 9.745721139753424, 9.224667932129877, 6.046467988579701, 6.803038333107534, 6.999878984619056, 6.785222830433441, 6.49364514784744, 6.72466944637663, 6.198476920581644, 4.857046365071446, 4.824458786254075, 3.8522204131427005, 5.066260699186752, 3.792330578055272, 4.409929679859648, 4.9225010266453575, 5.109061591582494, 4.566420157040571, 4.872431533398135, 5.26507415104937, 5.048051784170006, 0.29856362778010614, 0.7159308904623688, 0.5615184527370649, 0.0]),
array([5.81994281549717, 5.540673799382079, 5.183955333020539, 3.658358791088587, 4.361159065692311, 4.480037141192969, 4.002514622552013, 3.102625739835081, 2.991978674654974, 3.4930795603388956, 2.4891127196937624, 2.234067753661452, 1.506858355658958]),
array([6.531268195663432, 7.202783506550433, 3.7716091171243393, 3.686631741323999, 4.137168448979789, 3.666498409569374, 4.547368895698984, 3.895291093532664, 4.600422630893983, 5.11825866527672, 3.854064549460784, 4.489823074376847, 4.005040335374147, 3.8559529473446075, 3.296772126107151, 3.1617473437481323, 0.7415346576516213, 0.0]),
array([9.80374602717565, 9.148999994239539, 10.259209827744897]),
array([9.270321877126854, 7.89770197769618, 9.625139146308138, 9.065977591737681, 6.199188545326134, 6.387993867547573, 6.072003594694584, 3.635386769407682, 5.045126896175138, 3.8034159623561834, 3.464095794130089, 1.8806577623926266, 0.9205668432616533, 0.0]),
array([9.522799914034973, 9.311326942491945, 8.631725372407612, 6.030845911743124, 4.9838167727533325, 4.428412401672693, 3.8315112340099375]),
array([5.975201070757652, 5.803965844500189, 5.012682698875793, 4.360532022582877, 0.6779240312151609, 0.0]),
array([9.020902448025819, 7.445177505568866, 6.786940607601171, 5.6767014512943135, 4.263831075101135, 5.3006512800670675, 5.092653657965573, 3.8509844719057282, 5.148702270623459, 4.880535504086353, 4.5239755559316155, 4.92568358012295, 4.626996025598359, 4.011287841837758, 3.1716222733351795, 1.9031798789777072, 0.5693222756078435, 0.03849539265509842, 0.0]),
array([8.319339575848527, 9.216666949029005, 5.490546629200795, 6.029034381211039, 3.9112217769363182, 3.919103305963976, 4.6624970980423965, 5.167998574989512, 3.9733975596207607, 4.9558480410997845, 2.696630175961797, 1.9380103699326994, 0.0]),
array([5.604468232442686, 6.668852138814313, 4.380814074855867, 4.166700952521741, 2.8156044491139194, 1.698360175716815, 0.5260462526054692, 0.0]),
array([8.453516425637268]),
array([7.5147070869069985, 5.35888462341454, 6.917264302320347, 6.499516028261492, 2.941529816929895, 3.2905952439571533, 0.5303791205822528, 0.6460808718773501, 0.0]),
array([8.10225033145193, 5.419890375110345, 7.0190045704586685, 6.488208470746827, 7.128161576865426, 5.005134552958714, 4.801587009293343, 4.457661618380983, 5.226789980234692, 4.654263081920846, 4.443845474212312, 3.3990653527560153, 0.0]),
array([7.505313871393125]),
array([9.060622425727432]),
array([8.812899435329108]),
array([8.351923814791897, 7.314099794774073, 7.979472151035072, 8.57719712350573, 5.971973489572054, 6.061593834883225]),
array([8.047882500421942]),
array([5.8801083011356186, 6.1617052295821235, 6.612181582892, 6.920546445707614, 7.194856172248669, 4.071692103166238, 5.148307969686371, 4.817106483192111, 4.835758454339287, 3.638461722061608, 4.3777765038509795, 3.658373348697178, 3.124300608915891, 1.2758070674119937, 0.8915381826015232, 1.5133292217945522, 0.37401829693271466, 0.48836764533128024, 0.364958362160483, 0.44994419694647425]),
array([5.449893870373138, 6.718664407959327, 4.721608154576672, 4.369910354464323, 3.095267210939365, 2.0344400957855084, 1.858650589200277, 2.398976645547327, 0.46749847950533613]),
array([4.723628714400844, 4.103401983658704, 3.318538256355788, 2.713241089595617, 2.475034994519184, 0.0]),
array([7.563229343338847, 7.156862419156942, 6.672146720927486, 6.9162607921645405, 6.526630977899967, 7.0668957879608625, 5.963596724339402, 7.208266200545652, 6.6127699876317285, 6.088908572793781, 4.077221513287198, 4.056712657060647, 5.025011200720718, 4.165866087023557, 4.742513776774161, 5.141064824406074, 5.058591872493643, 3.5354614798653596, 2.9909383874266893, 2.8513654983230117, 3.473160400993666, 3.0484129459125566, 1.1278344549098345, 0.3377155811776425, 0.5500795945409577, 0.3212843724176511, 0.44933576795151264, 0.38370075583831575, 0.3754687012751228, 0.41649228788492676, 0.2591804152532484, 0.6718411539310895, 0.0]),
array([7.556789258046297, 5.578121439248692, 5.801648591419172, 7.228384177890802, 5.328648260676754]),
array([7.180722562280911, 5.394708248969687, 5.821239598030742, 7.193340101202566, 4.1563935004016335, 5.01096120894914, 4.3168556891976415, 3.632735821819259, 4.17962416641597, 4.517652141814137, 4.855805457985267, 3.9378294402212197, 2.835709339296187, 3.5129815133209052, 3.486251318258213, 2.402314406689763, 0.5045232481283846, 0.40925941251692277, 0.1585387758984681, 0.0]),
array([7.571809705459685, 7.338832270180652, 6.000854537976675, 3.651067031400368, 3.8710069557641598, 4.02200681665543, 4.007384578062478, 3.4796344652768076, 1.8336272247122896, 1.7747530572041004, 0.2739374767189896, 0.47830659433416156, 0.22294727955121252, 0.0]),
array([7.3785687803357805, 7.400574380267959, 5.52220656951463, 6.677357170620853, 4.545322537202078, 5.219552005530983, 3.661492392361378, 4.660051368362845, 4.515394927682911, 3.765077934025832, 4.841592551837345, 4.106703887309764, 4.673986145540107, 5.242857916229291, 0.4444923385584921, 0.13141146370279266, 0.0]),
array([5.689065218360325, 6.349571172381368, 3.8014509854418232, 3.0664227123590906, 0.521400014977761, 0.3327390942478352, 0.0]),
array([5.103401473689935, 4.777448203265481, 0.0]),
array([4.323010998763666, 3.940591038352173, 4.886419356006898, 5.022487406232692, 2.63022754077738, 2.9444564830536484]),
array([5.7315995589552875, 5.622363598933216, 6.451851812510725, 4.809331115438855, 3.773392957634063, 5.285031587722679, 3.7419976253734255, 3.924300564677874, 3.729020720046246, 2.7316172738433218, 2.908685996524199, 0.0]),
array([6.081661020129035, 6.209104319693345, 6.412523072632974, 6.068499051729104, 6.888887088583704, 5.566841521192904, 4.271387271664281, 4.3985831505759805, 3.8872317919359327, 4.473752579696728, 4.239489619782939, 4.025409766821266, 5.142139072313594, 4.957991162493657, 5.034269892052042, 4.760012654304915, 3.1637956814896064]),
array([6.46108475991439, 6.179573371527455, 4.4063573847735205, 4.059308638890847, 4.674702738472349, 4.861088888469327, 4.5203145148014645, 4.834904238042066, 2.95966200894296, 3.191267180964524, 3.4044299394865574, 2.355731735725843, 2.1202339343467282, 2.441060496941839, 1.595282217074509, 0.9406422593846417, 0.2943616419128819, 0.0]),
array([6.7680136890739275, 6.243854825133518, 3.7675177482597704, 5.275426087775375, 4.147665399428071, 4.785923236129212, 4.233618804025024, 4.444252274043688, 3.3139664969403486, 3.4690809149530546, 3.031039137722503, 0.3938972594392185, 0.6339967159989661, 0.7042972007618847, 0.0]),
array([5.475047054846373, 5.434362714370961, 6.424042619531867, 6.7165957765825395, 6.126279390490711, 4.962362438199794, 4.381733001345199, 4.094435929571333, 5.308809455078344, 3.399818762405505, 3.290211874023551, 2.116922280116841, 2.1536673582781956, 2.115407171545534, 0.49652110841412345, 0.6044683777869978, 0.41442891760804323, 0.0]),
array([6.203277057116599, 6.493999753657335, 5.384937148208595, 6.387221441760579, 5.406654779712276, 5.885837887083499, 6.243252325368142, 6.4481890947180185, 5.26484240100845, 4.856729376118235, 4.992969300084178, 4.9797444181556605, 5.215659563480437, 5.042185998194586, 5.266764749993878]),
array([5.5540474042022066, 3.778009916567629, 4.46818586410147, 5.05164745665918, 4.8459399501278, 4.239931448889974, 4.056294192622517, 4.848956715828664, 4.092880804216014, 5.232137403682511, 4.669207486139901, 5.141414112035665, 4.9442601654576634, 4.789059721629347, 3.5433315755175356, 2.9461925952281853, 0.7421919047295324, 0.6974024066964272, 0.40371838965926427, 0.4140968418933306, 0.3182860637248975, 0.39887853189859, 0.4616829206411607, 0.0]),
array([4.376985140174263, 4.593851893077808, 4.450745454482622, 4.6094127288828926, 3.654197800395788, 4.605644594845888, 4.438639497728343, 2.779489420798273, 2.0910221759248113, 1.1292050269099347, 0.8422179507737606, 0.25649123739044644, 0.5152589248696786, 0.4673503045831939, 0.0]),
array([5.385146641440307, 5.8828700126952285, 3.7221398993487633, 3.619914284904764, 4.671090869373333, 5.132594324674787, 5.152092119299201, 4.997705665517095, 3.5162476333200514]),
array([5.77827911013113, 5.6068883145947]),
array([5.65328453903578, 5.0860607556025, 4.37450002381492, 4.026122865848803, 3.1001224700588237, 2.4232222293420485]),
array([5.673255337510626, 5.6571288018436, 3.731695162954008, 4.196563706253572, 4.797634825957642, 3.9973123462436186, 5.172142897973389, 5.186976600456748, 4.95591158747345, 4.127234057113144, 3.750879620332814, 4.520058330177198, 4.386813817109705, 4.641529383243812, 4.622313000103426, 4.524148946054275, 5.085344321985679, 4.088876993895434, 5.144361668165558, 4.525514657290569, 2.8372399157102364, 3.1724253536005103, 2.7000907352928154, 1.0582542024880106, 1.5506726102482236, 1.0228894454388824, 0.6405628206711101]),
array([5.154482184044784, 5.213994533295605]),
array([4.670296885727586, 5.131291841486811, 4.306252452959243, 4.694404976373295, 3.8680133452208114, 2.822753161518946]),
array([4.994481740241291, 4.8004564320562215, 4.006560085746695, 4.839990019444226, 4.547185525463586, 3.422285012954426, 2.596553986266745, 3.2449293167626694, 2.1868751887834956, 1.3928100100619718, 0.4702757239016704, 0.7471204334056758, 0.0]),
array([4.298387811569802, 4.501152019838217, 0.0]),
array([3.9494829383039676, 4.212022686009703, 5.031860762306665, 3.966649203479517, 3.722979159014276, 4.812131666141388, 4.8467388388296655, 4.396451878357914, 4.511865193680964, 4.785463946776409, 4.496057669466116, 4.588932284505042, 4.288347001025105, 2.936826213988203, 2.1554282319631235, 1.894955686934471, 2.3742747904338124, 0.28933051957584344, 0.5032528711587074, 0.5299179316870011, 0.776825463255882, 0.44854132112293066, 0.5440748287069299, 0.547115752493637, 0.0]),
array([4.276837035785425, 4.370712569566714, 3.85988244153161, 4.96950713886227, 4.090789136750652, 5.012296518043422, 4.8714048584239, 4.8281700572869815, 4.335222425851047, 4.3780775704295545, 3.557852897874443, 2.646529595559587, 2.249523014938239, 1.103027316287125, 1.0956372817639881, 1.0220414095147756, 1.6920503739459214, 0.0]),
array([3.999279319704202, 4.506810165947028, 4.8324418453630384, 4.316503064539811, 5.0545575745643205, 4.688211252257937, 0.3885495960014368, 0.7286526191583482, 0.0]),
array([4.394859796827057, 4.395700306603187, 4.396308282636835, 4.734349579006871, 4.249844140462091, 4.60023085570237, 4.363103774727394, 3.634363493019726, 4.915327986746796, 3.8167326347660326, 4.1115258903848115, 4.476235195994696, 4.219239732538239, 4.3813537379992145, 4.914248153170588, 2.765987068404872, 3.0176523973315015, 3.5397722517178516, 1.9112674993391212, 1.2124497277394135, 1.2346235333771847, 0.6994792418507457, 0.5735394951407268, 0.44766341270026094, 0.738739341999441, 0.435919602582003, 0.0]),
array([4.9120257562694425, 4.2794819550436465, 4.23897574984405, 4.683958047477596, 4.13629219823689, 4.186294460788918]),
array([4.732877308186566, 4.20063211189995, 4.596088097686286, 3.8825003456520557, 3.124178299538784, 2.183268314963901, 2.46597791514337, 2.51645748634455, 0.8258193639524435, 0.2341769652729454, 0.0]),
array([3.60592800826586, 4.065999295315728, 3.8779543644680614, 4.626504004956901, 3.8260618459070805, 4.448349600820309, 4.160405151881408, 4.382088869676957, 4.5159974624024075, 3.318993417609248, 2.594173529253962, 2.2732058043480334, 0.7498081086502787, 0.5610892861553135]),
array([4.0560268150873515, 3.63379197232303, 4.271517038319232, 1.8442662042481301, 1.4536960362969564, 0.37859226154642195, 0.13422444897987773, 0.0]),
array([4.233009114164423, 3.716167163026544, 4.278556577377368, 3.6305937618328006, 2.805903007158944, 2.7894748151529756, 3.2326646141017066, 3.045313104709316, 2.7988522268671825, 2.344020160123993, 2.385329014253023]),
array([3.6721789447479516, 4.1420595687962525, 3.6030885164161566, 3.9665825587236174, 3.6442296838128643, 3.969560918808249, 3.0466858620057753, 3.012169410022378, 1.8183496364143297, 0.7688667432439585, 0.5032309329647693, 0.0]),
array([2.101517738921764, 0.40457919875588033]),
array([3.878822116918669]),
array([3.607869878904439, 3.495337364644267, 3.488539778059426, 2.6328179114512538, 2.3506484149993687, 0.7413465040102175, 0.019129174710607683, 0.0]),
array([0.8332093600570356, 0.16452186176904315, 0.2768854711389931, 0.0]),
array([3.7805668211230317, 3.850939778613279, 3.3508503829424514, 3.3987627756133523, 1.8684064453076898, 0.4232879259686437, 0.0]),
array([3.6643610862058296, 1.8561469492595735, 2.3931386489173, 0.0]),
array([1.841920049190469, 0.0]),
array([3.112229400062211, 0.15124813563260353, 0.14566847678861883, 0.0]),
array([3.2220324159023312, 0.5836139836881985, 0.2444558463297697, 0.0]),
array([3.3412174200799396, 2.4517516319593016, 2.015238081186717, 1.79534626657233, 1.537701093909096, 0.5716139110107733, 0.0216520593873491, 0.0]),
array([3.235936886189058, 3.1302386462832374, 2.511954012759753, 2.423031375361133, 0.44033871919538864, 0.3331038811908353, 0.10300043253872625, 0.0]),
array([2.806238314970959, 1.928081686758385, 1.9809871801536536, 0.6554455770596119, 0.6922795357139045, 0.0]),
array([0.1769708169443709, 0.0]),
array([2.3422895819644722, 2.3719049384635142, 0.14863905234188235, 0.6724547905952933, 0.0]),
array([2.725953689626627, 1.6167035773129017]),
array([2.6136986695275146, 2.944551338029464, 2.6006195364053175, 2.369711129898339, 1.9754869929751715, 2.221969268217399, 2.572158908436595, 2.2036372638455717, 2.4393400749671614, 0.8475397665723791, 0.6969075668788205, 0.13736318920422563, 0.17376513083944045, 0.4061539998799747, 0.3368438021758556, 0.27409317231457886, 0.3428776356816543, 0.0]),
array([0.8029496155689868, 0.6952127388599437, 0.7151276560713921, 0.0]),
array([2.8289768148652623, 2.370730556551054, 0.6282693720860933, 0.5419451120467462, 0.7142228025554136, 0.0]),
array([0.654392473402359, 0.16264409729051965, 0.6702756038989577, 0.0]),
array([2.605948571820911, 0.7702468589277353, 0.34620711275230337, 0.5399753051199571, 0.5042150349750041, 0.25410534915852956, 0.19325822548513927, 0.08653413118676603, 0.0]),
array([0.3062299634178498, 0.41339462262160515, 0.0]),
array([2.3452443469247304]),
array([2.282397859594835, 1.3234040247033354, 0.2401349487993819, 0.2867584638545885, 0.0]),
array([1.9097892518853155, 0.31758997692771934, 0.0]),
array([1.7217741345287787, 0.0]),
array([0.8400259889773527, 1.0698527414653474, 0.21653062761687458, 0.4406954567669038, 0.5060980732066453, 0.006501142917959127, 0.0]),
array([1.0948073334789763, 1.2104172376073823, 0.7889648575340928, 1.4049741044637059, 0.5970583528316591, 0.45551215607007617, 0.6186980488554319, 0.35572720549744813, 0.40888635507300514, 0.6762527992683022, 0.0]),
array([1.3566927670985842]),
array([0.9370784878644105, 0.2972696520506112, 0.5750019297015453, 0.0]),
array([0.2040520044680212, 0.6974446882527838, 0.23584531792979835, 0.0]),
array([1.057071083912505, 1.200462853481658, 0.3399688160816509, 0.0]),
array([0.8102577126231033, 0.32331744740055973, 0.0]),
array([1.1993404662362506, 0.7669161867240823, 0.43992144523517257, 0.0]),
array([0.7780748977700543, 0.692505036903705, 0.0]),
array([0.3264875525516918, 0.1871047318003809, 0.5173699238292465, 0.6737772399454622, 0.0]),
array([0.7558520922121496, 0.0]),
array([0.18948773598564062, 0.14722474057403911, 0.0]),
array([0.17981665709973949, 0.0]),
array([0.658489375830055, 0.2106256957535887, 0.0]),
array([0.19992584248854495, 0.3488730468912705, 0.0]),
array([0.2485516028749472, 0.3597848863742196, 0.26983440321236984, 0.24113026615740235, 0.5425009035735725, 0.2320887247388883, 0.19972731707750524, 0.12719725479109723, 0.16774544413250425, 0.0]),
array([0.4684967782847339, 0.5660078854833891, 0.6636080764597517, 0.0]),
array([0.132494409389798, 0.07464425519250134, 0.0]),
array([0.23639872883997504, 0.40715198238616, 0.0]),
array([0.39028296645173355, 0.3985353367155777, 0.0]),
array([0.4442618265714854, 0.0]),
array([0.23535979290072906, 0.15659543824383504, 0.0]),
array([0.2410346674988402, 0.4373622167593485, 0.20659351057156825, 0.0]),
array([0.4197448840014859, 0.41967049000001716, 0.0]),
array([0.13243645382124086, 0.0])
]
d = [data_1]
names = ["56"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T13', 'T14', 'T15', 'T16', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T29', 'T30', 'T31', 'T33', 'T35', 'T36', 'T37', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T46', 'T47', 'T48', 'T50', 'T51', 'T52', 'T53', 'T55', 'T56', 'T57', 'T59', 'T60', 'T61', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T82', 'T83', 'T84', 'T86', 'T87', 'T89', 'T90', 'T91', 'T92', 'T93', 'T94', 'T95', 'T97', 'T98', 'T99', 'T100', 'T101', 'T102', 'T103', 'T107', 'T108', 'T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'T117', 'T118', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T129', 'T130', 'T131', 'T132', 'T133', 'T134', 'T135', 'T136', 'T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T158', 'T159', 'T160', 'T161', 'T162', 'T164', 'T166', 'T167', 'T169', 'T170', 'T171', 'T172', 'T173', 'T178', 'T179', 'T184', 'T185', 'T186', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T195', 'T196', 'T198', 'T201', 'T203', 'T204', 'T206', 'T207', 'T208', 'T209', 'T211', 'T212', 'T213', 'T214']
def get_taxa_names(): return taxa_names