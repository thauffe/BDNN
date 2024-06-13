#!/usr/bin/env python
from numpy import *
data_1 = [
array([34.340795626645665, 33.25851426786267]),
array([33.01525196400821]),
array([31.35386488509456, 31.250689618378168, 28.035822339702108, 27.936032526843583]),
array([24.91821372558896, 23.93587136606404, 25.884290584656213, 27.13311791903342, 20.105631687960045]),
array([31.34585778264306, 28.457312728774212, 29.3492587748471, 31.43758191713025, 29.780065263115926, 27.604720296606523, 23.756402869325015, 26.014976620831014, 28.03847100881995, 27.678153067669857, 27.85011904453447, 26.867867058964997, 26.753469740248608, 25.793552406436113, 23.642353489606453, 24.968458978075827, 25.15857904535231, 25.901199799304578, 24.20921967118578, 25.107541218902867, 25.918689829490702, 27.10656848802182, 23.439529772471552, 25.762083166542354, 27.72702054206753, 23.217353280245927, 27.70631134131645, 24.6706901692708, 25.950983173793915, 26.860923193347848, 25.656630422441896, 27.386503362474993, 25.320910841468805, 27.09503505257556, 23.498864414055863, 23.88565105385866, 27.56145039409636, 24.805344745394194]),
array([26.46660230046061, 25.009132707818353, 27.40173829142021, 25.66578280871954, 27.816183934537595, 25.872479160117216, 26.947252116729917, 27.22928981921666, 25.59605440940598, 25.192880683652817, 25.91061614202126, 25.422804954550998, 27.615414091578018, 28.05751219770482, 27.17862467906487, 26.635197600228725, 25.92806008544031, 27.986332667429036]),
array([29.31004402673552, 24.102539375624495, 23.32216058309499, 26.6996514031481, 24.726252287513876, 23.635709750525564, 24.037858657982177, 24.084205177009018, 25.33627305687861, 24.79561000197815, 25.48518716108282, 27.698351401936556, 27.49312729831978, 24.442896640033087, 24.27667910946346, 24.0602739715835, 24.1250861672803, 25.76228361356646, 23.644215739823373, 26.88039327188968, 27.071058134395965, 24.549044849210517, 27.660188230183245, 24.715651111210857, 27.19000600234827, 23.833386108243964, 27.825655262534713, 24.07058174342163, 25.04267821910392, 23.946983047182783, 21.721574932227714, 21.69079871406376, 20.613396904771264, 20.442677715905305, 21.908614844546335, 20.539023165523425, 19.183458420711336, 18.31847803506429, 18.651062758405512, 19.151308730081706, 16.793509612612006, 16.297010527399557, 14.835987751133082, 15.853653758554602, 15.666412594805355, 15.154383822015955, 14.957721580651551]),
array([26.875773654956195, 24.14789910819424, 27.629858796138958, 23.665484901547025, 22.023800012769918]),
array([27.4169664503282, 27.609376970250118, 27.084910595146344, 27.34650204824255, 27.568852660858553]),
array([23.859742096784537, 26.277284564595348, 26.883386608673923, 23.26328959137129, 25.14680443584678, 26.72080692266652, 25.930627581247347]),
array([26.81618907045906, 26.447098763785654, 26.56452676421248]),
array([24.061646890857055, 23.92163140830028, 24.57755059443444, 23.39285642171735, 22.457948063360725, 19.203331761838022, 20.01922446301221, 17.83805818000416]),
array([23.958636320447503, 23.61791206923301, 24.065618224007675, 24.1982961217208, 22.527537054625224, 22.982457957815203]),
array([23.848746176230176, 24.102829872125085]),
array([23.442149476450854, 23.3966591133478, 21.622086124036883, 20.734603799894323]),
array([22.707432536259468, 22.20449786589038, 20.890003598883183, 22.05889458028882, 20.680769813581843, 20.584981407276704, 22.923462508107153, 22.645377657353873, 21.08662182607041, 19.808512155729876, 19.63521568301886, 19.65724488953286, 19.848723127197697]),
array([21.887177703652927, 20.611041193618288, 19.143575346157572, 16.7610246050478, 17.324354199425972, 16.410127826283226, 19.379599958208377, 20.153218577721542, 15.365143837906924, 15.88807957169306, 15.585730405933106, 15.80428834749236]),
array([21.359556832059603, 20.450844496409584, 19.320891771230105, 17.375153883302406, 14.45416491947137, 14.249722153090444, 14.796616544097958, 15.20119725949326, 14.842457435979767, 15.084492762925462]),
array([16.25239748822968, 19.54010519928116, 18.07042067295308, 16.14376381888537, 19.97595790411516, 15.315196923856922]),
array([20.648204938066282, 16.771539005408265, 12.080684108478541, 12.356617926350758, 10.688943828534871]),
array([20.546703047792235, 20.832602180232133, 18.803428702376532, 17.224287633864968, 19.674233372801318, 16.30498976845642, 17.888521247345658, 18.863809570602218, 19.913981724784126, 16.478934990244554, 17.785076547398482, 14.092688524590391, 14.391496086964585, 14.526599881066852, 15.466782489547084, 14.680633527874171, 14.143860773357762, 14.11728399664473, 14.803137333546186, 15.471013098295181, 13.958644510601012, 15.035424425586829, 14.48632259483348, 15.40114557093333, 14.458550304115304, 14.022521069309942, 14.98929786859856, 15.358611079363971, 15.554850524047486, 12.651416550779595, 12.489274662120536, 11.866031766496716, 11.881891143982067, 12.851698179892558, 11.787868731486926, 13.646374163523848, 9.893134476747317, 8.62028113494754, 11.492351633442308, 8.15897393736415, 9.315006947763042, 10.684884143403721, 9.382777164943931, 8.538502490662033, 10.218210319820091, 11.071448868048678, 10.085717758876049, 8.902073134955737, 11.430324978695202, 10.17733608150589, 9.199463400038539, 9.58818042214309, 10.950091601981493, 10.023237643415554, 10.39697550880931, 8.660742823910773, 8.473767391947852, 10.351029805819572, 9.373402305184971, 9.375911518878805, 9.140900425136216, 10.275535570217365, 10.781104288115824, 8.278238218316815, 8.841081534578743, 10.028576035956105]),
array([18.95456871413478, 15.89295034786974, 15.705833226834574, 13.117865811763755, 12.436977297194755, 9.425621950792927, 8.020035133508191, 11.053417846950753, 9.212392306867052, 9.972361992608668, 10.072223263257573, 7.256412442769648, 11.54281487038702]),
array([19.62237501824406, 15.27205510641403]),
array([20.532497494245604, 19.58618996296724, 20.37242414111024, 19.699920293991944]),
array([18.315326252972138, 17.42448515082846, 19.326908222335813, 18.72723473531891, 19.87470951419004, 17.36647991502737, 16.7941061447595, 15.176799094774847, 15.009774134078297, 14.192100836019467, 15.76904595955281, 15.795647932108327, 15.907340364717742, 15.344802744808643]),
array([8.475968500957176, 1.6675550252375586, 0.0]),
array([17.59533463439269, 14.16990347911465, 14.747805116403814, 14.318958559354254, 15.951922071159297, 14.548260457330281, 14.412620526997681, 11.95666791907003, 7.803281989020068, 9.872494202667086, 10.743034313565555, 11.236377510337103, 8.189335413920647, 8.328272199050502, 7.164519981037616]),
array([18.732923498043082, 18.536355158002547, 18.712078086766013]),
array([17.367913904824846, 16.453041214771208, 16.169992310976244, 16.48243423786781, 14.771742326929552, 15.373751169916277, 15.75394371127857, 14.712869611167546, 15.523566220166703, 13.886351183791358, 13.822695272022766, 14.965790218821812, 12.266293685378587, 11.884538712619715, 13.005076054041115, 12.780720820913693, 8.308892867654617, 10.981424015894033, 7.439247311134207, 10.543228642482294, 9.91154828902657, 7.29044090484005, 10.057018503163059, 7.5926986066132365, 7.413927994113957, 11.226688323703092, 8.226499504437099, 7.111389046391629, 5.3814262824304295, 4.80877196729181, 2.594198230200216, 1.5384492899646955, 0.6819787388375619, 0.07361738583070565, 0.0]),
array([16.11226319597414, 18.888172489580167, 15.55630046807078, 15.809070858261348, 15.559910995396534]),
array([18.994944574501456, 18.82515349588246, 18.584226216368144, 19.1302181272945]),
array([17.25875805863399, 17.501287365549913, 16.830705258998833, 18.75861827198822, 18.875606572395327, 18.76699582812342, 15.418175128012217, 14.14407832006638, 14.74793735281015, 14.727886670475385, 14.981010685416958, 14.630628134627665, 14.498234244306634, 14.646810857205642, 15.69125203562619, 13.907696636627781, 15.546432365784524, 12.351686547707688, 12.070042085614197]),
array([18.753253979617543]),
array([16.37513480539835, 18.36304445368465, 17.548142908318773, 15.296933958441869]),
array([17.539277632647103, 17.896297461258172, 18.519537452321384, 14.288458445912838, 15.598726224192545, 15.039124532088255, 15.00065461564343, 15.03827357589883, 15.533369508520337, 14.318149035932269, 14.97880317132964]),
array([16.327116506185916, 16.711243765804067, 18.23779276701475, 18.378372251158925, 16.821571234607738, 14.802147329009886, 14.00967088483135, 15.43187646140236, 14.332566986832951, 14.4249174870891, 14.830858326159202, 15.965423534223788, 15.890864069853878, 14.404843157617398, 14.974654776751525, 13.94454890276591, 14.906246222561714, 14.997580011696112, 13.748214616826289]),
array([17.287499343375245, 16.536555537600965, 14.626341196240517, 15.721731555431257, 15.179634879510855, 15.392365155529529, 14.667594097195746, 15.185966353928809, 15.206292526486708]),
array([16.56043975816888, 16.49372056315815, 15.129531548604756, 14.841030320764283, 15.212924505859231, 14.895748699816806, 14.914738886364113, 13.888050883139881, 15.394243945272201, 14.966231583460118, 14.48666864462183, 14.461482761748114, 12.127397813841032]),
array([16.197885791855615, 16.355453876354883, 14.868308157366096, 14.30682939143051, 14.46040975946377, 15.240528083578706, 14.905232027510273, 15.472910876768518, 14.246965405685641, 14.132880589561502, 15.815790080194635, 15.914716751507747, 14.436836423393707, 13.013054006396247, 12.211604302374901, 12.148900967716402, 12.313107172882471, 12.293647450561824, 10.438594247798694, 10.83556544075351, 11.173183502311323, 10.995812984787975, 11.584671750655488, 11.218783875334376, 10.875984793853064, 11.403984290425612]),
array([14.872616366360344]),
array([14.608401988913986, 14.346802063315934, 14.639586323972027]),
array([16.93806401816406, 13.961902569240952, 13.706916425009293, 12.255924892397653, 10.936969989773885, 11.526439514534886]),
array([16.7494193543089, 14.586170252951241, 14.985996765654521, 14.706639950460247, 14.450880684020866, 15.485278363942351, 15.281073097557123, 14.83469970965731]),
array([14.623226150941433, 15.041142417483279, 14.790336418382667, 14.686058931365231, 14.063805316741316, 12.834590995223934, 13.41220578318569]),
array([15.74670034746113, 15.7953746683485]),
array([15.971694988594736, 14.691866147865428, 14.923131646666963, 14.449398913012216, 15.356182385049864, 15.275328347167788, 14.782214146335479, 15.579960631254124, 15.877429184835123, 15.173933269641703, 15.105145401207901]),
array([16.034978673300895, 16.098065754343544, 16.060038403210093, 16.04935574108083, 13.88306226862172, 15.547313757445707, 14.537601812561965, 15.111740204666814, 13.858102973203351, 14.19277600268363, 15.17644715961651, 13.982063917909885, 15.375094389174949, 15.422355911941452, 15.030189998578797, 14.636451416553593, 14.857337937161141, 15.94527586006706, 13.905059340057809, 13.932844798308992, 15.382616833846095, 15.633243379014957, 15.155670568050395, 13.932322275037606, 13.927323020264083, 15.593855423494688, 13.961444885462242, 15.197359412885111, 15.4760132061409, 14.150719830790878, 14.327304880955579, 11.780085434287438, 12.030878863595076, 13.742524910533065, 12.68970788288076, 11.710268608986006, 13.027398360986243, 12.829889177902574, 11.349628780528718, 11.60594701482488, 10.24263885223987, 11.28663676296422, 10.206026910909841, 9.779381352337227, 11.113363875161754, 10.27124457526167, 10.850965034853829, 11.30793289273865, 10.32974501397703, 10.741501773534877, 11.402022603763058, 10.957938628083744, 10.022080349226622, 10.906305643371136, 11.471104212475629, 11.590091636762187, 11.333852924762358, 10.363839801384561, 10.965944883974693, 10.903382607629633, 10.902211997118172, 9.902214904491709]),
array([16.19283699275043, 15.36405528425576, 15.492828523275708, 15.831869932852703]),
array([13.87876936588755, 13.859166547418795, 14.63561076791575, 14.897826484502747, 15.051565064349468, 14.92183842253937, 15.782657875846796, 15.748883188714466, 15.287490144178907, 15.498069163424239, 15.140266681174207, 14.841623516027406, 14.204800978398724, 14.461578117007033, 12.31259998351172, 13.677141186270543, 12.482330123840883, 13.496313946032236, 12.263565898226636, 13.034749069282082, 13.60126137114333, 11.994026220733861, 10.463558105577267, 10.930172760881394, 11.31251854644143, 10.499257325794087, 10.353291978128425, 10.36809868671451, 11.231261790859216, 11.493608935810672, 11.08850929423722, 11.564886460158572, 10.861348014528653, 10.430739232257205, 10.328851013480778, 10.922101716684404, 10.358415947340776]),
array([15.306880862533362, 14.252096113212533, 11.904786849078855, 12.218254823831913, 9.8964835549395, 10.57630015365059, 8.600510333349797, 8.308664711981617, 8.38183827402081, 9.316656676807398, 9.343813340970552, 7.54542439946927, 8.656752652964016, 10.392696823032345, 7.520158938962771, 11.12144114531951, 7.675780947904462]),
array([15.08912325643605, 15.14711733088023, 14.266699795614091, 15.197442834256933, 12.481609633247944, 9.789964737957146, 10.239630381563655, 9.743539100110718, 7.506743185417442, 11.421108998398521, 7.837616380216754, 8.804141628659066, 7.258854253079385, 7.350239234736587, 8.726790667706808, 7.3742020133100326, 10.095987502208313, 11.246074503868103, 10.729554693455114, 11.41089046160833, 11.471670536745041, 7.942320264581416, 8.74844086247479, 8.63976967895962, 7.705195836892738, 6.2851450190637985, 6.335184040854235, 5.520295008104824, 6.474358218017266]),
array([14.005658850867555, 11.907211319046485, 12.95761177342643, 13.772591559009811, 9.702323399530842, 10.566188166426347]),
array([13.947137836575354, 14.108132003609425, 11.632396081330427, 12.08665828416716, 10.619212070633886, 11.334783116134915, 10.048109057311278, 8.63423630629031, 6.437776244669463, 4.92508990154527, 3.1758246328685638, 3.46163720429311]),
array([13.915756695583429]),
array([14.440653975123388, 12.534777703515314, 9.200701071342225, 8.239233063763354, 11.382299921869997, 8.4613708534217, 9.017703151416024, 9.83310309044988, 9.32767841106271, 11.229995548662524]),
array([14.028607861371391, 14.014306428302717, 14.28494648877368, 14.429932147815784, 14.252851140343957, 14.25659829220618, 14.222533289504517, 14.095772905828467, 13.877323741789596, 13.502165728480673, 12.222956681354887, 12.964237051278834, 13.14609141135072, 11.589711373973426]),
array([14.394358400651777, 12.836282457247224, 13.543290477329597, 9.785120687513999, 8.799634543677797, 7.6516319227547225, 7.285631377600032, 10.850798966767194, 8.097963768819262, 10.215815355046438, 11.5942044259125, 8.124583116781928, 7.792118035766368, 9.53069403042752, 11.168023344115756, 9.812044895839337, 8.251225260707734, 9.91296322644125, 7.120678111350168, 5.6880892495542845]),
array([14.044000819993936]),
array([14.12968964880243, 13.429723132380312]),
array([11.800613699308514, 13.213579728331515, 11.30555228165347, 11.612928641368804]),
array([11.02594642362022, 6.592776196688195, 4.1848311331533985, 0.0]),
array([13.616756246682115, 12.358770355460273, 8.773995121474737, 11.15250768788047, 7.867115918070137, 10.647282427508555, 9.098711123751542, 8.185138939593402, 7.621424777277135, 10.886909695017938, 7.696503550998425, 7.6569472058975645, 10.475623525244508, 7.116002500885223]),
array([12.797519818565226, 12.258875105725103, 10.441974327066967, 10.891044243296854, 10.951932512399077, 10.779135941877136, 10.503182672977886]),
array([11.81333060310447, 11.93018618340428]),
array([12.823127576534297, 12.335226753268103, 12.838048582358264]),
array([9.157066838159473, 9.139587822660161, 7.86854454794615, 7.371798212853119, 8.444617053099263, 10.621652322738397, 4.0112488101920425, 4.851993747109155, 3.0051332648592433, 2.001470110475542, 0.0]),
array([11.65268016437217, 7.822668818926335, 10.855620741017717, 9.927110261266474, 10.84990459712046, 11.557680327043117, 9.598065335449611, 8.861795215128925, 9.849066138304266, 9.696041347816156, 7.4729410825517535, 8.860316986459, 8.858430788661739, 8.92162337142863, 10.981962572338619, 10.154723710075887, 10.60065272907367]),
array([11.335464518322608]),
array([12.099106243900634, 12.49736887352497, 11.433204504755482, 11.051205059150947, 10.8399594461863, 10.972010854291506, 11.041686415447817, 10.493188050010351, 10.733285836172131, 11.299870084965324, 10.33859715483482, 11.071406044726846, 11.161921926283421]),
array([11.89350787419692, 12.240212621299863, 9.311395371271441, 11.172655331312633, 7.666574853267261, 8.955508079862167, 7.745183129182292]),
array([8.444844132214136, 11.618484126365098, 11.32395735224346, 7.833251305612597, 11.245069534881527, 11.500061961278469, 10.374040659548118, 9.903525285643084, 10.473644925023615, 9.453385126922313, 8.289460119555105, 7.2774521234083505, 10.828661334142076, 8.796665542592915]),
array([11.732887632496759, 10.583232219069124, 11.391489758337094, 10.870258729635278, 11.013804430821018]),
array([11.647507665117312, 10.462415798158181, 10.577672310653591, 11.113775845454605, 10.41165878018702, 11.0972913106584, 10.661336380211505, 11.219915938978861, 10.54713548693547, 11.462917882850197, 10.052858715914619, 10.89169023713544, 9.545141655109909, 10.168210580513701, 9.952809149016833, 10.09121590470923, 11.224718751227515]),
array([11.860295481994283, 11.024150644308103, 11.299759809882346]),
array([11.985534031568449, 12.082641047464673, 8.303467385504268, 11.346968225312196, 8.673989062869907, 10.185181829369274, 7.617622695756437, 10.742256007416401, 9.491992393455234, 9.889114346032905, 7.678636167585554, 9.779309274049645, 10.95036315221698, 11.152992072182714, 9.784721595656922, 9.341711519758995, 8.520838661080433, 10.034460540817467, 9.595712102616904, 7.520987578143128, 9.65620721159557, 10.788689753233099, 9.75183202945373, 10.898976145061217, 11.048045929492575, 7.528855944661346, 10.13564758460946, 8.042328980019029, 10.249307727272242, 9.775154688406543, 7.641151224861177, 9.819786170631076, 11.260841178467189, 9.917812459385651, 11.542697331135251, 10.355400836729938, 8.840404579183891, 7.686837820138756, 10.957458529557316, 11.554202009258773, 10.094095439937139, 7.419489349119045, 11.42106694950921, 7.870373835980472, 10.209437643069084, 8.824263395278344, 7.395610405996854, 7.447455354964576, 11.147799748571948, 6.272312814683182, 6.909078459249489, 6.848846364196887, 6.895878593917549, 5.619217001404392]),
array([8.981753901070853, 10.812155209748232]),
array([9.142512539338723, 11.383369866908918, 9.056822119276367, 11.540033340425293, 10.488231247023055, 5.70020255769276, 5.728706449202586, 4.819753689175656, 2.2153191656280584, 0.2950924811321591, 0.0]),
array([11.139804312612988, 8.66493002538811, 9.238418920005818, 9.410301963916647, 8.389998961961243, 8.11076209024688, 5.042870823166917, 0.0]),
array([7.7949702496938205, 9.192193011765623, 9.468014663038796]),
array([7.298493841754385, 8.528778508770875, 10.923381788666207, 8.630612147508348, 10.362627482951712, 8.045121603251724, 10.097810600556345, 9.983134472714962, 3.712110894678151, 5.010857030916031, 5.290034009527614, 3.3829902346862553, 3.469104211455467]),
array([10.798148953807567]),
array([11.10148664245919, 10.764992714402029, 10.374069305001974]),
array([11.088791460609468, 10.955686475894112, 10.472577408964614, 10.723035452001085, 10.604224944416647, 10.738562462134588, 10.453215845672913]),
array([10.216308488421934, 10.074179789939889]),
array([10.755088218018951]),
array([10.544317258265481, 8.158775360737092, 7.559219802188629, 10.260299052522697, 10.331239604660741, 8.889245518400191, 8.439313303141583, 7.3373118076663655, 9.545705235278934, 9.06952496405983, 9.356342287463658, 8.859025613239403, 8.722068500574178, 7.479321452643706, 10.545975047670339, 7.663905752500289, 9.728108314364379, 7.547877957587265, 9.893386166541966, 7.781803421455811, 7.578149669220061, 8.66515052116658, 7.595043453162701, 7.297617063825747, 7.779597083305514, 7.9278587497900235, 10.020306190477767, 10.144849761680517, 10.876713695644481, 9.26823382771943, 8.738204330027266, 10.26083740321812, 8.246431330696037, 7.667825842718281, 9.875966323130758, 8.157167838949668, 8.208576025466666, 10.745651420219776, 7.59076103283024, 10.777948772301471, 9.850846982996206, 7.392959305080623, 8.25979025668838, 9.008068237061082, 8.310275998205437, 8.163947986215714, 9.894943380460612, 7.725877864391853, 8.498963395468525, 9.436806693930622, 8.291475452723795, 8.561039463445358, 9.49896218366559, 9.686348138655068, 8.340357605783911, 5.772019178943948, 6.342348348265112, 5.6901808878561475, 6.554066998002638, 6.0672149007976195, 6.960254698369727, 6.740194970072534, 6.51152592751857, 4.049615601926882, 5.177217073984436, 4.989199227754111, 4.597339590041317, 4.810752592579599, 4.06715286519299, 4.484789722163154, 4.439621267185371, 4.81188951437012, 2.903117307325885, 3.5884785763390896, 3.2507400481150035, 2.9048688695767764, 3.2962053680367065, 2.439461673285585, 1.4783090825091085, 0.2156604883316271, 0.1529177643913503, 0.7150190941178551, 0.009204549311283117, 0.08842694563682384, 0.05492081851663684, 0.11721148051537754, 0.0]),
array([10.640503650229931, 10.569103294262517, 10.355487588901555, 8.872008834056876, 8.356577042637246, 9.244346269909816, 10.592427106827968, 10.400701488793922, 8.472953436472233, 7.670256896061499, 9.838294636632266, 9.990389587735866, 7.925571196532129, 7.853190894740295, 9.94143669973093, 8.457880277697207]),
array([8.789536845497803, 9.440670081041374, 9.15050466949251, 10.693459155944229, 10.544010101764574, 8.403664457250484, 8.805053990544204, 9.200815076594193, 9.48535184031875, 10.199046822880984, 10.232818819740356, 10.209052891491744, 7.901770406804406, 9.115529839945772, 7.520393096237973, 7.503367644131343, 9.004028443252489, 10.489740369647723, 7.736524215799644, 9.010892796934746, 7.311254097462549, 8.9427306519986, 10.154612004054805, 7.960615903171258, 9.889458759751363, 9.935807606267247, 10.570840678457818, 7.894142449250223, 6.075323421847122, 5.85290248852908, 5.6786366766866365, 5.973959781133424, 5.365329008154145, 6.964050400575379, 4.831527736060792, 4.016115447324349, 5.124321105259806, 4.474949168793554, 3.8168390029349957, 4.538280381772053, 3.903953761083745, 3.9095188855814924, 5.017191347234165, 5.046941748648907, 4.079812448100366, 3.677558840160973, 2.741558618467843, 3.479474004083171, 3.515917491857738, 2.812964935102772, 1.5179327853690383]),
array([9.97069925858084, 10.63579895504887, 7.359983387641042, 9.32566960235442, 7.412829200337297, 9.68476337991181, 9.320150711736272, 9.171016756206333, 4.752992446436368]),
array([9.760751143717268]),
array([9.527647391728275, 8.015750919172328, 8.784673023617929, 10.167062897892402, 7.9398296707848885, 8.279434355359186, 9.282695158931348, 8.23374012138002, 8.253802669312954, 7.739794465918435, 8.178108837524508]),
array([7.522586021050176, 9.051811978869907, 7.634559278825879, 7.851016062125958, 9.654465601671955, 7.788802216602324, 7.885258287481393, 8.435236811351356, 9.195695937711006, 6.243797065191033, 6.445517702003009, 6.080731634739336, 5.118489630743881, 1.8359946059402625, 0.0]),
array([8.014811732010566, 7.277570868521111, 7.558520924755335, 9.088088344818368, 9.437046119606208, 8.68525924175441, 7.834647950325076, 9.30820206313808, 9.353578993125627, 8.94083355340612, 8.899485826739076, 9.565551508038192, 9.119728951807005, 8.645038711351793, 9.63018364874837, 8.680642179645494, 5.896050079933295, 6.0490554755612, 5.91379959828759, 4.707547096530788, 4.476901996349308, 3.118100648026222, 2.267617399793095, 2.2469293624340354, 1.2600258534981408, 1.7229087532488412, 0.39776396398323904, 0.670629395604809, 0.0]),
array([8.978739713719033, 8.72706338041107, 8.51747535630381, 8.045084378326676, 7.318356456451196, 7.74221705811827, 8.841764297960307, 9.480071549852838]),
array([8.893407497283732, 7.727813854556948, 5.938972680633363, 6.516830200888119, 7.237720951618594, 7.106715591671608, 4.11806665303425, 4.028499187071056, 4.7040948514814325, 3.5839560946582676, 3.107660127371155]),
array([7.564125347157105, 8.744694457213432, 8.879424284519985, 8.855237111575219, 7.539455002045646, 7.584463866005784, 7.473601919197623, 7.838230033160507, 8.199320962880224, 8.39085737827179, 9.223338816892692, 9.06553874547136, 7.4715321735920135, 7.611991419738919, 9.641315466984748, 8.391201627261909, 9.00745936278601, 7.471044470602851, 8.070514548958544]),
array([9.233943763797937, 7.945311811806477, 7.454362328667021, 7.496410877599679, 8.696962452134365, 9.172617229697824, 8.454035624815575, 9.05770795585848, 8.915004422754484, 5.445059257099423, 7.12437450526748, 4.673448165380913, 2.3719535747737956]),
array([8.456977569152228, 7.611567651471494, 8.808725811389687, 7.382031406151386, 7.785967651948224, 7.544156299246902, 7.800638104721705, 8.824405429966017, 9.041325340048513, 8.875844234608602, 9.047921644116085, 8.624006293069153, 7.247980547841577, 8.658183037818134, 8.366570241179957]),
array([8.357622919872494, 7.878368635354068, 8.61527108711043, 8.581155391193404, 8.523174323653333, 7.865822535647206, 8.410267670840321, 8.077090484750123, 8.308696250826712, 8.38428244798601, 8.346269988990967, 5.970856762896549, 5.714683148501599, 5.952015574330551, 5.219149766491288, 5.166885355292886, 4.341403064719175, 3.4300550986183334, 3.2155928390228357, 2.9576834238215217, 1.2810895345728697, 0.9181863557592476, 0.0]),
array([8.993840152619196, 8.62041673977544]),
array([7.905154378320826]),
array([8.500621438368475, 8.564264720165637, 8.876347065866666, 8.012030927168807, 8.14477363237841, 7.344427403147407, 7.405373446776069, 7.4598789455273256, 7.581816693118835, 7.400710759691332, 7.412695439049235, 5.894180367505315, 6.053202427408733, 6.706374151515082, 5.946391227938156, 4.955545842329444, 4.092416535026253, 4.194309715463005]),
array([8.625665257270926, 7.492146770225769, 7.5273637879195165, 8.272418436398691, 7.427451876808595, 8.012642915535936, 7.929822006611685, 7.48078795002106, 7.726156624738534, 8.252160716566168, 6.306852069867784, 5.895767884879005, 6.462409305353376, 5.9371968021733945, 4.399527106527195, 5.314031003129462, 4.988413812738793, 4.193019632193899, 4.212621097375036, 3.860789280999068, 0.11829275497894412, 0.0]),
array([7.476595009098086, 7.404672977656575, 8.452662143340737, 8.184876450417242, 8.849610496655997, 7.242484161772656, 5.837297464909766]),
array([8.729156749841598]),
array([7.866653583901582]),
array([7.531185468144939]),
array([7.734945754863446, 7.406371958246866, 7.296855137588763, 7.10043460640255, 6.725109279350908, 6.872702838715145, 6.6161516923953005]),
array([7.3210514552367085]),
array([5.661524772398689]),
array([7.519149778156169, 6.86508544243258, 5.669352484980612, 7.22910295160579, 6.775757742500101, 6.652722201057191, 3.7935452590797336, 4.631011617553868, 5.062419419279118, 4.724715728595448, 5.046442791974103, 3.8715399239602606, 4.5304383718681445, 4.61382275903485, 3.704393217604424, 4.8664524323546505, 2.930118789024898, 1.874799177233601, 0.6047732281358631, 0.42715200561843814, 0.012981408236947226, 0.0]),
array([6.284121318836048, 6.518710927728684, 6.191706051421825, 7.037010657913048, 5.2436820838235025, 4.26336049205402, 5.059856095472039, 4.405655432254623, 4.4490171424488745, 4.95798000947209]),
array([7.3731930923462965, 7.37386859951494, 7.342712859536592, 6.837804858415932, 5.968347626085096, 6.672295813958335, 6.625884387675819, 5.4932617152441345, 6.182769923808628, 5.425309129720082, 5.050956193490098, 3.873059085760579, 4.983751516628632, 4.17001437882961, 5.126489270561982, 4.839727645926084, 4.704268199173762, 3.8155353491285138, 3.387851067936544, 2.781358295274226, 2.4213811983133406, 1.6420541504446287, 0.21527720989838206, 0.3370206246052709, 0.49659697499530725, 0.27989481937769933, 0.5210407420758593, 0.7472259945368186, 0.7180268209031314, 0.0]),
array([7.290351725744873, 6.126702475682975, 6.69725861605632, 6.861696489019096, 3.931372853726992, 3.688601303007278, 2.9951056751860827, 2.2786336186033807, 0.3806494661639173, 0.11675161297001065, 0.0013446631290800065, 0.048735981670199324, 0.0]),
array([6.231075193641619, 6.720127478961852, 6.117434421318522, 6.541047680833916, 6.0357917363484574, 7.063606135581628, 6.414767400429225, 7.073544329379428, 7.135155368216601, 6.372375752397176]),
array([6.342607032637803, 6.360112036087951, 6.5811027277397, 6.817719904486938, 6.550729396015066, 5.324522077693121, 4.751884951237102, 4.2328048684373885, 4.44300318651628, 3.773007311733141, 4.93065317665553, 4.854659233280171, 3.2873598702910773]),
array([5.660807187046359, 5.3485566805619404, 6.375607121387921, 5.085432785331671, 3.9951979946864737, 5.144695725056822, 2.7056348223154307, 2.001371794061428, 1.3573731324970024, 1.603350331880478, 0.0]),
array([6.551183211418278]),
array([6.13091427455851, 6.250003880526798, 5.614210943935331, 4.315789280689952, 5.120143205419151, 5.191808112846773, 5.119530026815056, 4.806278624733551, 2.813229316517631, 2.8999686698791827, 3.4862649805489068, 2.6933356666218677, 2.6869290784964175, 2.0678162277937524]),
array([6.626690202718619, 6.453185395597774]),
array([5.027877625300314]),
array([3.0193042647712476, 1.6693372621506173, 0.3494093841616684, 0.0]),
array([5.867396606523306, 4.928295418770244]),
array([3.988353027409108, 3.6658382305996415, 4.320303796564401, 1.9982906540887204, 0.6775715972688472, 0.0]),
array([4.427560791589828, 5.1187002246524616, 4.3992857631927285, 4.037753842403722, 0.49230631672012887]),
array([5.692082461393566, 5.024714735263359, 4.8870534083962, 4.288856179512456, 2.7317602548154234, 3.546318736185668, 2.620621272656873, 1.8864710271821075, 1.224404925398506, 1.6465293855732241, 0.703573880297168]),
array([5.5295255811808754]),
array([5.707673305642652, 5.622982922085473, 5.584669558121384, 5.2090410194854515, 5.198484446301502]),
array([5.484163377941713, 5.390162530650266, 4.128209342592219, 4.865918529037892, 3.9770301161220867, 4.901487314963172, 4.143352752141283, 4.703218742901961, 0.5373972189540251]),
array([3.908073662146422, 4.379343694755359, 0.683459656373426, 0.0]),
array([5.044912280932359, 3.4322127366859796]),
array([5.4452596128749615, 4.436272437511877, 5.287787117589342, 4.2892614685219845, 5.088252785737074, 3.264393561926256, 3.538579148891118, 1.6247897599630305, 0.6608220614000997, 0.0]),
array([3.8993648456009287, 3.7897333304632195, 3.2951811625661733]),
array([4.595335327075202, 0.0]),
array([4.961059537533354]),
array([4.673480713798175, 4.008267085860038, 4.906170491819624, 2.875583282638102, 1.6399475994712185, 0.2302420061285908, 0.18565443626874223, 0.7226820473968569, 0.01121587295675762, 0.0]),
array([3.728547755150635, 4.542320902857191, 3.9013496587373577]),
array([3.8219568231964303, 2.8039567177030302, 3.552864298088757, 2.381404328901169]),
array([3.767279578819771, 3.3109702673037305]),
array([3.6648815252996814, 3.6258143551196067, 2.894346564166477, 3.1152666306919308, 3.417379969694009, 3.0771310815892754, 2.9445560469716794, 2.3261751349857063, 2.5569846597470702, 1.8209918914559857, 1.3610595948381987]),
array([2.844223471482385, 0.0]),
array([2.6852437593370464]),
array([2.7953241437682532, 3.463446163058598, 3.0762635497906032, 3.309919617116349, 1.2247619137680832, 0.3715529826980704, 0.13594865707405834, 0.0]),
array([2.697427951911357, 0.024926236365861845, 0.0]),
array([3.4279504176299143, 0.4657976420340483, 0.0]),
array([3.242783574386845, 0.0]),
array([3.3346101548864286]),
array([3.195649746158267, 2.6518117731933915, 1.508600236891615, 0.887836905162916, 0.45671086136533534, 0.05127112543334379, 0.0]),
array([2.775728655178942, 2.5832324951382213, 1.0558065869312236, 1.23929237480629, 0.0]),
array([2.935564774382644, 0.0]),
array([0.1682060625143612, 0.0]),
array([2.577390452620168, 0.8723891336019625, 1.4803685755105132, 1.2511404335447884, 0.4601690089936936, 0.04045563482115938, 0.0]),
array([2.140637503492068, 1.6685421694379736, 0.39545696146778414, 0.24441516473798208]),
array([0.40425087302932633, 0.4096202672727404, 0.0]),
array([1.7922962071158266, 1.2998506942401875, 1.4919180393189249, 1.0314669612695981, 1.7369314720126285, 0.7255328797967127, 0.50129911149868, 0.08268546628699024, 0.0]),
array([1.1178663562258233, 0.878345858660757, 0.8237962662057906, 1.0601168624184787, 0.5607789085853394, 0.5971267252483301, 0.0]),
array([1.2962785974694921, 0.0]),
array([0.8334911160791124, 0.5212287841787557, 0.0]),
array([1.480304952395113]),
array([0.6797406648721934, 0.06315665844336095, 0.0]),
array([1.8518044318989653, 0.7949667502230382, 0.010869026746032698, 0.0]),
array([0.8171010376174519, 0.0]),
array([0.3160137283791902, 0.35289179753834615, 0.7742384413935722, 0.0]),
array([0.32924738753393457, 0.0]),
array([0.5551848931833943, 0.0]),
array([0.09128767128294019, 0.0]),
array([0.14939647951321833, 0.0]),
array([0.17628651774540627, 0.10265223863014813, 0.0]),
array([0.18315891924222605, 0.3529947612088543, 0.03184977092284537, 0.0]),
array([0.16142719325058635, 0.0]),
array([0.1934703562003112]),
array([0.038355917631811516, 0.0])
]
d = [data_1]
names = ["48"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T3', 'T4', 'T5', 'T7', 'T9', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T23', 'T24', 'T25', 'T26', 'T27', 'T29', 'T30', 'T31', 'T32', 'T36', 'T38', 'T39', 'T40', 'T41', 'T42', 'T44', 'T45', 'T46', 'T47', 'T48', 'T49', 'T51', 'T52', 'T53', 'T54', 'T56', 'T58', 'T59', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T70', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T84', 'T85', 'T88', 'T89', 'T90', 'T91', 'T92', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T100', 'T102', 'T103', 'T104', 'T105', 'T106', 'T108', 'T109', 'T110', 'T111', 'T112', 'T113', 'T114', 'T115', 'T116', 'T117', 'T119', 'T120', 'T121', 'T123', 'T125', 'T126', 'T127', 'T130', 'T131', 'T132', 'T133', 'T135', 'T136', 'T137', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T144', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T154', 'T155', 'T156', 'T157', 'T159', 'T162', 'T163', 'T164', 'T165', 'T166', 'T169', 'T170', 'T171', 'T172', 'T173', 'T174', 'T175', 'T176', 'T179', 'T181', 'T182', 'T184', 'T185', 'T188', 'T189', 'T190', 'T192', 'T193', 'T194', 'T198', 'T199', 'T202', 'T203', 'T204', 'T205', 'T206', 'T207', 'T210', 'T211', 'T213', 'T214', 'T215', 'T220', 'T221', 'T225']
def get_taxa_names(): return taxa_names