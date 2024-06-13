#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.39073008696652]),
array([30.458975891308548, 29.707963757468228, 28.291294022189664, 32.56360185819336, 32.06669863591692]),
array([31.87348777513336]),
array([29.429824650566562, 31.202229206197458, 27.186583715156278]),
array([30.11996324792123, 30.0602530872902, 28.215734324826336, 30.418743212728106, 26.950442345934675, 27.42491407488612]),
array([29.45820827455411, 29.116615483313844, 30.17152865880373, 29.17994355822275, 29.485765419563805, 28.735110521817294, 26.457424105607657, 26.902750067564448, 25.76490899400854]),
array([29.002614909169896, 29.798161507503877, 27.644748923923217, 25.438159410863232, 26.654277424511243]),
array([26.38734602519676]),
array([27.77225030652608, 28.084329330819664, 27.96604930140001]),
array([26.555357242264623, 26.74903512504714, 26.60847621799856]),
array([27.48799000321402]),
array([27.16322442509387, 26.28245216328058, 26.092922699744783, 26.23514384884352]),
array([23.920989884540653, 24.419613358859223, 26.562323183014648, 24.770248770309664, 26.420082825819705, 23.690920030085636, 24.105373617039895, 21.980535581190594, 21.275332563843257, 17.623650140123022, 18.08466779709473, 18.186716883020253, 16.97809962050157, 16.471668502236152, 19.967995184952773, 17.016680920642223]),
array([26.541636666744775, 24.62212823542052, 23.27066767809596, 24.64846097952369, 21.610992689057724]),
array([24.700859028037844, 24.819873525597405]),
array([26.09581546778204]),
array([25.3085825527015]),
array([25.37485249219946, 25.114584439382078, 25.71890300957553, 22.37947168345525]),
array([22.24257578226243, 21.64853132512365]),
array([25.696068528966872, 25.1615547947013, 25.12228219003249, 25.63564529634651, 25.481139349264517, 24.204160172844155, 24.806934371223605, 22.376630914407105, 21.917806930147453, 22.876370371495103, 21.310537957549034, 22.2418067352571, 20.66748163427569, 20.182430274228462]),
array([24.16324164892125, 22.645157818115536]),
array([22.29982719337855, 22.079067205200115]),
array([17.169086823346255, 20.217480849299683, 12.376185557209883]),
array([21.416351397391114, 21.99401223333385, 22.04836363894307, 21.76581347089736, 22.409774970687852, 21.97830786849963]),
array([20.87055416858371]),
array([20.7909633204211, 20.72633080283945, 20.54054052865832, 21.533538804764355, 21.74507208524453, 19.811398170820905, 19.16740180859476, 19.981933469085543, 20.025129742696745]),
array([20.98970055182115, 19.77171193734928, 18.397219089265608, 18.488748840800838]),
array([20.2164691317907]),
array([20.91868611487901, 19.88722908218226]),
array([11.624529268836504]),
array([18.072904950661233, 19.152969744932147]),
array([16.103305156756246, 18.283965473553724, 18.913427682838602, 17.145595749539382, 18.18938327202733, 18.827421499311193, 18.717499188824974, 16.119541377828114, 17.52435712984729, 16.672091586125234, 19.172593940845495, 15.915504640160515, 15.408998623607545, 15.52518461255748]),
array([18.891035948148527, 18.194422678722464, 18.10262464505572]),
array([19.203169251121192, 19.257819470255427, 18.091689818609446]),
array([18.41904043415475]),
array([17.951760181170968, 17.99873386982895]),
array([16.15090974221495, 16.62741260070622, 14.82224976329281, 13.778798337169967]),
array([12.545471657303413, 8.490038886859574, 9.233009790834824, 11.531029734049497, 9.362666029790283, 6.258715682659895, 6.267860003373269, 0.0]),
array([16.204929356273308, 16.950109486164592, 13.865417817446346, 15.39195657145317, 12.661269049646034, 13.436843932133682, 10.8502947642339, 7.555760069229132, 11.264681116932197, 8.750172653181139, 9.541751957234922, 9.580398327845723, 7.9716419324266745, 11.122373715540395, 9.646553894508818, 11.045886817004783, 9.346668326632381, 7.150173665710876, 5.502610529115119, 6.295107028950575, 6.3937377952017425, 6.5912428935985385, 6.7709083331853455, 5.905550121543713, 6.0298390374419455, 6.4078497651622675, 4.214797754914431, 3.807613254144801, 5.158012602901518, 4.144626828826247]),
array([17.188668042894236]),
array([15.28167638726034, 15.266586080258914, 15.458854751919423]),
array([15.637202062108482, 15.525772046294115, 15.924159111977394, 15.195736547799985]),
array([16.20489156235683, 16.102890332123234, 15.229049071142644, 15.050911372006109, 15.804013509855809, 15.290791491651936, 15.514045356611106, 15.268240211506878, 14.497608029742157, 15.884013167587913]),
array([16.40212599352033, 16.301104977474722, 16.520510220454213, 14.321824018092578, 14.270985994194206, 14.411209604931239, 15.364020227473903, 15.069947739519458, 15.744749521435303, 14.888885091938482, 15.107484734524999, 15.057593514855823]),
array([14.27212688592537, 15.10518476419674]),
array([16.23770157115185, 15.280307948580383, 15.36541912433979, 15.3602918734406, 13.821132248300888, 14.901151187467256, 15.191179155281059, 15.315205976051816, 13.87708156402244, 14.10375479194773, 14.969736354046821, 15.426882455446721, 14.02041345435593, 13.62997167070997, 11.860874006201554, 13.369579466604254, 12.54895413349349, 12.581715039999539, 12.689474415243906, 13.6505353348402, 13.31716789300584, 13.139078151075676, 12.957876953469098, 11.314104506402154, 11.079868605498026, 11.537262042260084]),
array([16.324603138334428, 16.294293467084092, 14.317012001874273, 14.830268599662647, 15.021156680283262, 15.02483332953762, 14.407234385659033, 15.904881894061928, 13.92993119720542, 15.8686360861584, 12.503777504463777, 13.000020383877066, 12.331347185338002]),
array([15.300268736010722]),
array([13.863651439278273, 14.150014714574588, 14.255485114258839, 13.175901915674192, 13.465562201554654]),
array([14.821500194714327, 14.757281739714445, 14.632496892654451]),
array([15.510549485004137, 15.033578091915688]),
array([15.035678889400257, 15.009748152521013, 14.63014438626712, 14.140715973490007, 14.165858672385566, 14.089720607739025, 14.150097133884977, 14.131590028109292, 14.568914148701605, 13.119576612895749, 13.707372075365633, 13.001676863737863, 13.819704777696646]),
array([14.711113489311142, 13.661433144157838]),
array([14.819956085375049, 14.78617998567077, 13.849942158444453, 14.473304319803027, 12.264319170494346, 13.28374729641823, 12.034880913834872, 12.081849556504235, 13.636560694741174, 12.627602414225898, 12.716663224133416, 13.728224474925984, 12.593254573636392, 12.654635972769517, 13.598591473922497, 13.58255745167923, 12.777541641795668, 12.680466266757985, 11.340089545958346]),
array([14.227276803720185, 11.64496830605179]),
array([14.166201897899237]),
array([11.075586281369484]),
array([14.112582415602928, 13.471686474736245, 12.984402767236727, 10.906833633717195, 9.400372681986976]),
array([13.948966573027304, 13.588886269112123, 13.645445059820924]),
array([12.287448031874526, 11.981944611951262]),
array([14.117008586762797, 12.05479248717181, 12.151887633806483, 13.502739845038873, 13.12172626518993, 12.244366211762246, 11.900161765976295, 13.34804123840646, 12.466576444833244, 13.018797278011878, 13.434180501017998, 11.365999233631765]),
array([14.033120038316492, 14.052827182211876, 14.046585244918413]),
array([10.159733749724348, 10.820178825700255, 10.595766006218858, 10.260908900904592]),
array([13.928101899985277, 13.880737378451284]),
array([13.758307954308753, 11.760297130681819, 13.34786359356817, 12.557946908393571, 12.7637539489632, 10.890105757882683, 10.178852240037296, 9.508014740110951, 11.297434870642586, 8.302969831564479, 10.511945047425803, 9.209505566235677, 10.731429533864757, 7.506609855459388, 11.065056398027668, 9.925509553254129, 8.166474063806398, 10.27039339144768, 7.478528733409887, 10.859535253298418, 7.725888270257041, 8.530193382068164, 8.469334669362304, 10.266852636984392, 9.943062607355168, 11.192439424596758, 9.768378989886763, 7.5899090902266835, 6.7386875371640125, 6.917027914267849, 5.6002234613726145, 6.536273242655443, 6.047154486295655, 7.091634994288755, 7.024761808916458, 5.698357011484229, 6.294362520485117, 7.1723041315134575, 5.663405472292796, 5.87189045672151, 6.828997292997615, 7.135738121353905, 4.44049889854619, 4.482027670400798, 4.684971790182807, 3.344885010204452, 3.5524365003800056, 2.897507555227654, 2.0139034130802314, 1.705912204167506, 1.3189419241013327, 0.9754989809451438, 1.2081042679732628, 1.0217202276887476, 0.6885716670160076, 0.2702853119015872, 0.0]),
array([13.303594205766444]),
array([13.55365711851192]),
array([12.964926405080352, 13.131219039065545, 11.63564405033589, 12.363870791278801, 11.936312783379464, 11.873217495999906, 13.554893146045021, 11.574554515997656, 11.413594929188452, 11.288201147810813, 9.02483656594492, 9.42950721532324, 9.146257644687738, 10.482204108825092, 9.896561794690111, 11.018948220993217, 10.642594628646558, 9.275095494636222, 11.369517131624264, 9.525186547131424, 11.03843022644686, 10.189293432157882, 9.355987620648534, 9.652256565802757, 11.161641707051464, 9.127311647469412]),
array([13.25905245045844, 11.68640348746252, 13.679111988168021, 12.034417036178723, 11.920492631128973, 13.378907897655234, 12.670083621511319]),
array([11.712803396155643, 12.01508458099727, 11.951941895607604, 11.900294681948058, 11.684090290633863, 12.024461741604437, 12.123605245156995, 13.082669882975715, 11.375931642156168, 10.755136196784738, 10.159399327942207, 10.149988298535348, 10.284073289990843, 10.355408821522985]),
array([13.371702434760905]),
array([13.371034778424695, 8.842973849239192, 7.381509311256093, 9.603585408004026, 7.452615157849285, 7.380993068661517, 7.485726816180768, 7.539644406707452, 7.141034119012927]),
array([12.825716480543456, 10.89311528284565]),
array([11.805281697515744, 11.823565623769769, 12.410479577816567, 11.687009589877794, 9.618333862581638, 9.80870718306759, 10.927862691461932, 9.916996276278201, 10.254879265370437, 7.331586891113124, 11.028251484895478, 9.526986193524086, 9.49123959639878, 7.583926606921637, 8.769078988977874, 7.753043162408834, 8.621761449824984, 11.316788715181106, 7.503483644276376, 6.812295935460712, 5.478240046025832, 7.208102053139223, 6.616387002028307, 6.646497999413612, 7.195716770329909, 7.18559715918848, 7.10877067519536, 5.917498129913324, 6.2814565186158795, 6.4592261411340415, 6.948528386595238, 6.225190269668242, 7.120992669002647, 6.599643431351517, 4.538150488452646, 5.302774026519013, 3.8906128764062284, 3.947886631613188, 3.8389740949867255, 4.998612285329335, 3.03031375774915, 3.3026592508269403, 2.8074585830818792]),
array([11.300915307833161, 11.447897574595512]),
array([8.786998506162035]),
array([9.391743143734317, 9.054951645904545, 9.007349504259585, 7.514902973030338, 9.62321355558481, 9.28598362522257, 11.246204457499694, 10.505542164740861, 7.540710698758402, 10.163290326734332, 10.487422023068586]),
array([10.60286647092397, 11.340102381877141]),
array([10.41580909973164, 10.204581132890493, 9.568683058729466, 9.715156992668295, 10.397178430203493]),
array([10.591093672293225, 7.6138166467046045, 7.997051511380731, 8.44039573605688, 7.857989680955487, 8.793723216571697, 10.79108019913406, 8.573326024567478, 9.646679594086713, 10.604243561563877, 8.872697956359932, 8.402511754722127, 7.876780824422734, 9.296474136010696, 9.967346899003367]),
array([8.070375361336545, 7.359220793839116, 8.432719490073834, 10.638265463397676, 9.27993641350874, 9.461648381305265, 9.474916672862383, 8.091596687672087, 6.6297965235133764, 5.232054758722701, 4.638504725229623, 2.610615994423531, 1.6598583534538067, 0.9030448447823378, 0.0]),
array([7.643081952756459, 7.353809544669922, 9.432552593394929, 10.003588458425979, 8.278959622906736, 8.098525222925376, 8.803164908280522, 7.469818173694254, 5.622205147068759, 6.386788588753478, 4.73914230265672, 4.797595100800622, 2.74742756885863, 1.1122701130810808, 1.0518679712293775, 1.6530440543727796, 1.3220074000947393, 1.367565175238119, 0.0]),
array([9.34389571508328]),
array([7.275441459770465, 9.32427076547065, 7.933392473203506, 9.306569927598515, 7.835388296988926]),
array([8.000101880328069, 8.584889490268395, 9.785659665884683, 8.245680695198685, 9.529710299408128, 7.36004842539212, 9.503346457006023, 6.480288166693998, 7.052638449638745, 5.568602246170606, 5.918349396599717, 6.6527378598501965, 5.776748478632896, 7.203537095070844, 4.65753606458008, 4.371434332549198, 5.213260623209107]),
array([9.8496966814674, 8.188916385041228, 8.057079756181441, 9.987638242427883]),
array([8.447674397013103, 8.228353406102816, 9.053597433844576, 8.36957478038045, 7.4637980630281735, 7.907373581887578, 8.082802427704376, 8.801591452365672, 8.485475671053832, 7.799106342511673, 9.671287741096972, 8.970117321196952, 7.332616717974656, 7.319856571273849, 5.80637188031618, 5.553496499443121, 5.608235272526063, 6.549487323469004, 5.4382813581715155, 5.8135596541790875, 6.272770421044554, 5.740001759861905, 6.566950757867064, 4.783206936899532, 5.252851267679043, 4.8229997080022775, 4.834529471800822]),
array([9.5085072223047, 9.75473410590646, 9.620324384623286, 9.751583896548917, 9.705720237252844]),
array([9.17454142381208, 9.679445641060514, 9.05797537255159, 8.242721942508652, 8.951014605488114, 7.731677158515979, 8.091352968393117, 9.490073463681474, 8.988643119156105, 8.14280245138357, 9.428760835693687, 8.978977494644937, 8.87445858389355, 8.831191372730064, 8.388353215242017, 8.16485562206693, 7.9379772312350525, 9.47928069124965, 8.187612823634293, 8.263169754740789, 8.249857364609667, 8.455508034094999, 8.713456037968282, 9.107727779630052, 7.814685003288819, 9.098625820773517, 9.061351610047822, 9.52532289487444, 7.770234351023775, 8.151324487903064, 9.372711533906443, 8.698725612033762, 8.013348842796262, 8.653521764593414, 9.681669080655439, 8.562641695839973, 8.695167760484283]),
array([7.949697625496767, 8.670311313482681, 7.394642953955504, 7.7086273662075335, 9.389325209159152, 7.7458857167126425, 9.501130864015078, 9.434691739900755, 8.44951027393936, 5.564105665398067, 6.555383893451575, 5.5132683374544325, 7.216741640276689, 6.062679493681026, 6.412396640195869, 6.059450384649349, 5.769831694602048, 7.139961291121181, 6.903392080945729, 6.820934228879985, 6.218868199477026, 4.6565522093497735, 4.561938808745005, 4.126301181972304, 4.784201659902806, 3.4763198384128455, 1.8499357160987164]),
array([7.510377450827614, 8.331421651151484, 6.240164901248251, 7.219877599306882, 5.4105859537254135, 7.176065824722526, 5.731500494676284, 6.208286686923123]),
array([8.339622714833993, 9.47882061621374, 7.334972540930589, 7.542148832981184, 8.679668964208055, 7.763694787777613, 7.429920974640976, 7.760556820848157, 7.980273280101674, 8.809443128541277, 7.831024161196007, 7.908427194248806, 9.022783366784154, 8.69389603392952, 7.51682391407563, 8.787228221675967, 9.034979445928428, 9.229331538065312, 8.50817129104151, 9.549879771132309, 8.65480388819359, 9.274618180754679, 7.6877129443941525, 9.135029531577098, 8.418441416130907, 7.69905220748132, 9.321763869550344, 7.567324256974069, 8.908352019706728, 9.24244253375545, 7.127705166633691, 5.473669333979194, 6.7650500828841755, 5.887270328651061, 5.953600093619734, 6.227958331718311, 5.843985665740404, 5.5362942874352985, 6.262135878467162, 5.604789272519228, 6.344287998479204, 6.593604622010512, 5.963138383095735, 6.119501718374249, 5.857396510950739, 6.682955828397075, 7.057384890673649, 6.929922982503637, 6.958420474993954, 5.702086912636028, 5.65933515393532, 5.412329461898582, 6.511568760191322, 5.43193380081126, 6.59021951630414, 5.064897395325972, 4.137750588979351, 5.056151976742297, 3.9649825897014104, 3.9056740117188546, 4.543756215114236, 4.609631695603073, 3.708215239623396, 3.2106321110900424, 3.1654862224501805, 0.8932742550625837, 1.1422794422727418, 1.6090291198535498, 1.7360975669153733, 1.4716356056372444, 1.2102830381795955, 1.5322266910484164, 1.5056426353206307, 1.5154732282107934, 1.0199647827340117, 1.0697471872797666, 1.6571553878480545, 1.542463614134235, 1.3552626258344787, 1.5500616991034115, 0.9267606083937142, 0.6654019591407444, 0.32956631075297854, 0.049662942114302736, 0.03132108106618729, 0.11373879000453921, 0.0]),
array([8.806927439619548, 4.5767500452924885]),
array([9.202814908515476, 7.485069410678408, 7.3257666102849175, 9.271973114158346, 7.445916832006464, 7.321211124596482, 8.848157636152171, 8.873459828430668, 8.704764522806256, 8.095719899749856, 8.309238568804039, 8.910028634845261, 7.587704029874285, 5.833624129311377, 5.875006820892704, 6.24156741411381, 6.708997646528623, 7.0005126053440145, 6.060671793114583, 6.453885597094882, 6.299261740740929, 6.578799592671753, 6.537060998871548, 6.175111192282716, 6.322610161489603, 4.123996932502899]),
array([8.195166147763262, 7.197063110107535, 6.91339905207545, 6.622670998038347, 6.763508659855391, 6.8264443708311555, 6.925072151218514]),
array([8.841004999192384, 8.224110678054075, 9.180531269768117, 9.296694086998933, 8.609243560524241, 8.490788014832415, 7.6858523098698175, 8.84867498913735, 8.348012884357914, 7.768282053791736, 8.461896081732027, 8.760090653540766, 7.77520116499374, 7.324445548525527, 9.348561695412139, 9.104292286564919, 6.819961818557755, 5.802518778853497, 5.529105423258547, 6.931463574971306, 5.924597215029748, 6.458358430578907, 7.214175436428954, 6.911050172423985, 6.18980829239024, 5.898497020754128, 6.85575143547669, 6.176093839751853, 5.755612004391841, 6.444530315059569, 6.4017119207116435, 6.510882570229622, 5.573543784086686, 6.968761832028254, 7.103179731404682, 5.676114904546099, 6.834775680254298, 5.133214501114617]),
array([7.895565289020257, 8.703957455716807, 8.902890764472316, 8.904131794924995, 8.609835967521132, 8.417057220587301, 7.939111544461893]),
array([8.034807125839995, 8.39312783085937, 8.911986288683348, 8.9520354642904, 8.51255168282851, 7.366833927413852, 7.755362430135637, 8.256611863734811, 8.331030620900085, 6.07214948359124, 6.652017370314495, 6.479354940988089, 6.598339367503508, 6.972843814608656, 5.805675823982518, 6.388327592773195, 6.422922530447776, 6.628571921426046, 6.755406739194008, 6.026102044181267, 7.101288759118446]),
array([7.815839663869371, 7.894322820108583, 7.3248320275115315, 6.843546559381843, 6.2484473418862, 7.037884492326214, 6.749018454399919]),
array([8.094986443542568, 8.168469107087224, 8.782295316024197, 8.420301959088142, 8.121928723567693, 8.16802504140938]),
array([5.933180856313564, 0.0]),
array([7.972423400215696, 5.706122414087356, 5.797576436899353, 6.927914248034752]),
array([7.633419114052256, 7.5863312274250685, 8.539307798871645, 6.444960872399637, 5.531109827593216, 5.202346253907473, 4.792672467711017]),
array([7.341601786538609, 8.70954757379384, 7.316987773396171, 8.459630374606135, 8.305398029027419, 7.262774425413359, 8.316778469235855, 8.064161712264497, 6.4573368091531504, 6.858979216891949, 7.053984213844573, 5.783544644123924, 5.953894393625394, 6.453862104975055, 6.634411328327113, 6.804078020045182, 7.118715927398876, 7.063106869651547, 6.374044293369096, 7.0030159679078485, 6.774663831267451, 6.304097226461476, 7.199915622534128]),
array([8.664638564949087]),
array([8.518185634212152]),
array([5.642756909255491, 1.384577973771868, 0.49894509161957584, 0.10693742477565854, 0.0]),
array([5.943563784302275, 6.22079705771842]),
array([8.135559915681936, 8.433160475195237, 8.27124561433194, 8.393608758904467]),
array([8.28130604757217, 7.540408906319547, 7.778896647535046, 7.276769862503903, 7.325203479306326, 7.796127251818827, 7.532187555596119, 7.594416916197435, 6.537858828563527, 6.526267932703819, 6.985898378168693, 6.926668526121016, 7.0019341073768055]),
array([8.171297783622517]),
array([7.638011364336178, 7.5929201329519795, 7.507381968077869, 7.532519215522429, 7.94254145380153, 7.261799567861276, 7.181841275617549]),
array([8.026368379651144]),
array([7.468129646728353, 7.522990158788254, 7.325842614071041, 7.611427035855735, 8.033695034438352, 8.117719497614528, 8.066989696941453, 6.7109033420255235, 6.89431654543511, 7.078938148953466, 6.053200312855619, 6.3513004482006945, 5.718468217945603, 6.072352036355102, 7.052937558405934, 5.828041861030679, 7.194715710680345, 6.207028924576853, 5.581843419933022, 5.60125676687563, 5.679516473597281]),
array([7.270636086698939, 7.358976553504803, 6.505505122861906]),
array([7.848027033873688, 7.293254729242933, 7.641353652179924, 7.249335138920507, 7.895297173986338, 7.75331237134073, 7.7703485823526535, 6.485670629380628, 5.963802738430363, 6.052862395509411, 6.44733594360029, 7.053834347591461, 6.051212848755079, 6.894553023271539, 5.951415617117331, 7.052292385950277, 6.113324746862691, 6.270055973189665, 6.872234465940404, 6.689409373682794, 7.233423497455506, 6.081162844966569]),
array([6.81198899526825]),
array([7.624126735354493, 7.678828327767824, 7.860987847113679, 7.493101569217303, 7.275204958630585, 7.443120960335414, 6.816791262866859, 7.183097876809887, 7.215624008657794, 6.930894055549457, 7.175469978995506]),
array([7.260968308639835, 6.995466651643969, 3.2656787257413793, 2.7223691841141515, 1.7197614438316446, 0.8648926582524684, 0.6900253221322259, 0.0]),
array([7.345688106827014, 7.435769942563645, 7.507285097841013, 7.508073180046126, 6.173608572846906, 6.848119080966555, 7.049553161099795, 6.539932477425422, 7.0821991554691, 6.909925720948989, 6.665238376083715, 6.167243992037718, 6.835995084014213, 7.088983566785985, 6.711955539204125, 6.118745758521577, 6.581981004455136, 7.212703967533752, 7.020126299208378, 7.1037826435027425]),
array([6.930929189003942, 6.210643902618885, 5.8187034646149085]),
array([7.297963739315449, 5.998018149921661, 6.250828993011214, 6.4506065731593365, 6.774223809830007, 6.947761492285396, 6.636108835163124, 6.798042519975171, 6.95925083126978, 6.248615467882934, 6.064789748298306, 6.401799054213293, 7.151668344317078, 6.282419868222476, 6.302970194470568, 6.850384863131307, 6.493084713233449, 6.307133286843122]),
array([7.033273523643657]),
array([5.460660062803254, 6.284872043716156, 4.283052820819478, 4.9300096244462335, 3.6915003916136784, 3.4548782929360833, 1.028644496406415, 1.5419025760362335, 0.0]),
array([7.105263944647837, 7.065342651807526]),
array([6.814806297682447]),
array([5.589910915163593, 6.463909284517827, 6.773938403880896, 6.87470959539067, 6.248091471621578, 5.430573485245655, 5.4906174789954445, 6.794812012653725, 5.459105778283426, 6.166405023367202, 6.619570384533559, 5.700682602795805, 5.738968361385856, 5.461383288759282, 6.479018721115126, 6.733961603637335, 6.544914037818208, 4.36082468644282, 4.10317403172421, 4.057247328712587, 3.1684385198776712]),
array([5.6538381794964705, 6.67879061678844, 4.858865901952565, 0.0]),
array([6.146915945075149, 5.4791407614864385, 6.36344467175648, 5.569618889966165, 6.161563069269189, 5.7603537951584425, 5.854919788183527, 5.485104173937898, 5.6419181684587985, 6.000944840832927, 5.448710099763808, 5.9236273149398215, 6.2879149516824, 5.824599982085624, 6.432801451154035, 4.59728485023009, 3.896375989349317, 4.967793093018881, 4.773868113588442, 3.810359494031145]),
array([6.662690378714337, 6.297317750962898, 6.104657327644777]),
array([5.3820064230979945, 6.544221451816117, 6.712322552935887, 5.727987481663654, 5.724836408618629, 5.5023823866613615, 6.187890785470922, 6.398859049714197, 5.527288527772353, 6.174391355898889, 5.380370412665698, 6.625703617113092, 6.049302887109436, 6.148982243974093, 5.311229360757078]),
array([5.470443353412967, 6.335293494882092, 6.396061345998205, 5.85735845241298, 5.68933481278631, 6.4653008628079585, 5.727542018112876, 6.518872716398867, 5.3418057725512424, 5.384363703576171, 5.644435794703415, 5.4595212381747125, 6.5700168041198115, 5.103108624301995, 4.46882652376562, 3.736568513598649, 3.8991016346588196, 4.721432469857783]),
array([6.581333821902822, 6.4854071362333805]),
array([5.619613485908575, 5.587447428511271, 5.713086229352519, 5.921402579805644, 5.916570413602516]),
array([6.190450915682084, 5.959359748115331, 4.974242361321117, 4.553283195916724, 4.674489729368668, 4.0382420050559364, 2.1331559540919045, 0.0]),
array([5.957065471726617]),
array([5.70911216487566, 5.443015393786121, 4.144072115413341, 3.192500468717589, 1.7524366687844204, 0.0]),
array([5.842354251681665]),
array([5.410943710408176, 5.023875608112233, 5.0681107187084535]),
array([5.938849943544253, 0.0]),
array([5.932847687959834]),
array([5.647553353873963, 4.823544922470677, 0.0]),
array([4.55866193919476]),
array([5.716876172680359, 5.468337455408987, 5.431675112608401, 5.740959564362658, 5.844500548279506, 5.33222818914891, 5.338382882488041, 5.530717707848817, 5.67894614591746, 5.559590258365393, 5.512447917482098, 4.585122291041426, 5.16373214742746, 3.78792289201793, 5.0420883079973695, 3.7535101600051917, 3.8543952310899416, 3.8262482596550207, 4.17259310750568, 4.359814845455265, 3.656968013323641, 4.373784470661957, 4.86491998091329, 2.6568834247702364, 3.1845041076810254, 3.082220424478289, 2.9047598077297625, 3.2761669384349426, 1.808190455068727, 2.2082580944613106, 1.7941228228441473, 1.6773971898071303, 1.5711242629174267, 1.6912243031430814, 1.6839254925661167, 1.6454505214640984, 1.6996447554738876]),
array([5.3982122189498565, 5.027256056432985, 4.750014382973508, 3.074764816259706, 2.890400132617162, 2.5248284963833125, 1.71773214605474, 1.7650944594544382, 1.408719334033095, 0.0]),
array([4.672061236969813, 3.1606773945955]),
array([5.407950811700669, 5.025278139887361, 5.19261471462179, 5.164512020063704, 3.4232866415454257, 3.2998699623359458, 3.30862872913003]),
array([4.231572848421786, 4.936361036262985, 4.882206434761512, 4.184640521841932, 4.952370445506681, 4.901263144527354]),
array([3.8673914015473, 3.6989601352886106, 4.990047960990908, 4.430877581519197, 5.06946604042683, 3.6981326257078098, 4.945762159396405, 3.853652867002054, 4.919312351907886, 4.6384715478542216, 4.347819478067439, 3.9624822508784177, 2.915576547059385, 2.94050653524812, 3.4933016336696543, 3.288548859529646, 2.9028417648890685, 2.40701113159244, 2.5555016146280845, 0.9816580979165567, 0.9105364223247486, 1.349307066591764, 1.2108647246754067, 1.6185310331113836, 0.7995802833522005, 1.2380300988342006, 1.0842309442501594, 1.6084648112013236, 1.2255133612255973, 1.3699931334941287, 1.2284921194177865, 1.2196719253441946, 0.7796611592683754]),
array([4.487157132010956]),
array([4.230415684582462]),
array([4.211624667459063, 4.294975534368919, 3.068699628621208, 2.44481287958606, 1.3806119006622046, 1.6475569398150929, 1.6741488310587485, 1.5486365955064851, 0.27378748276128406, 0.0]),
array([2.9225526362642738, 1.120290265682145, 1.4783443605428972, 1.2482394415518214, 0.6656194550713748, 0.0]),
array([3.5279505075699396, 1.4755221381093806, 1.6957615230368759, 0.0]),
array([3.8804131343874664]),
array([1.4614585270189693, 0.7990293720440655, 0.0]),
array([4.015013255000647]),
array([1.7155082381411995, 0.5884204747395076, 0.0]),
array([3.371768071536569, 3.3803770411690053, 2.5903158457681714, 2.9228744556971606, 2.209357784354483, 1.6797658098726223, 1.3999222232587156, 1.1544533862368094, 1.6384172006631086, 1.7544004402126023, 1.7414010508133746, 1.7912393881948885]),
array([2.207715002072198, 1.406636546167917, 1.756994137704017, 1.2697783917998144, 1.750208163876561, 0.9746002802014762]),
array([1.9129561934169321, 2.3165357183621405, 0.9114018292802435, 0.9292154068626497, 1.635312802403427, 0.0]),
array([0.09080728329164331, 0.0]),
array([2.8474273817696063, 2.5724056674813687]),
array([2.0357354444360283]),
array([1.9536976335499139, 1.4462896751083745, 1.3936868200437558, 1.6205389286882734, 1.7903357119799799]),
array([0.8823189486387852, 1.4485327471865674, 0.0]),
array([2.6566675126670947, 1.7801537167776023, 1.2939018017303388, 1.1629605750982996, 1.7697417491532157, 1.5742838396204377, 1.5781671974736713, 1.5635465265986617, 1.2855983272580909, 1.3990394584670733, 1.1318998883067501, 0.44270342508283606, 0.6324061236831724, 0.0]),
array([1.8229952446615245, 2.027504066250931, 1.490566616088864, 1.5762808511010893, 1.0686615036726321, 1.5082640228758795, 1.692283445951113, 1.6930208485614235, 0.930250011505801, 1.2183391716752734, 0.9552781324915784, 1.7798016725388803, 1.0220362323806365, 1.3407806188777518, 1.0011127131153499, 1.7690873471945856, 0.35383551981594896, 0.0]),
array([1.4575618079120012, 0.5825569575259582, 0.0]),
array([1.8128918963815182, 1.9122374699492433]),
array([1.0144143682687667, 0.9680463190584376, 1.2562398911358696]),
array([0.8630943947474861, 1.1400040665627176, 1.2070576812078815, 1.6471294248019064, 1.334200096626911, 0.3264054056288475, 0.18097172384096538, 0.38095892217338406, 0.27355278283324846, 0.037493238137682544, 0.07527306238053995, 0.0]),
array([1.0907804598249957, 1.0121419711767625, 1.2229293710086406, 0.0]),
array([0.9969648733637146, 0.0]),
array([1.0422132218881295, 1.2468577958733957, 0.8823648507315278, 1.0613569102128981, 1.126905554615235, 1.2850301642955912, 1.2381910715510147, 0.9203753211096696, 1.2151714417285358, 1.2136853185399428, 0.7321881577564024, 0.7046196759529681, 0.6733090590418238, 0.10755576572980027, 0.0]),
array([0.881166790591951, 0.9338677517747629, 0.8886090279033511, 0.312748673710281, 0.0]),
array([0.8090686634178783, 0.8059844296707711, 0.0]),
array([0.8331764804876022, 0.8098602203786083, 0.4848071602508811, 0.26852759487058997, 0.35368251958603547, 0.06725167857952646, 0.035906735205245385, 0.0]),
array([0.4459223763937051, 0.6321190414922604, 0.11938882764844334, 0.0]),
array([0.400567067993857, 0.0]),
array([0.2537810919748702, 0.0]),
array([0.029338811852864313, 0.0])
]
d = [data_1]
names = ["57"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T1', 'T5', 'T6', 'T9', 'T10', 'T11', 'T12', 'T13', 'T16', 'T17', 'T18', 'T20', 'T21', 'T22', 'T23', 'T24', 'T25', 'T26', 'T27', 'T28', 'T30', 'T31', 'T33', 'T36', 'T37', 'T39', 'T40', 'T41', 'T43', 'T44', 'T46', 'T47', 'T48', 'T50', 'T51', 'T54', 'T57', 'T60', 'T61', 'T62', 'T63', 'T64', 'T65', 'T66', 'T67', 'T68', 'T69', 'T71', 'T75', 'T76', 'T79', 'T81', 'T84', 'T85', 'T87', 'T88', 'T91', 'T92', 'T93', 'T94', 'T95', 'T96', 'T97', 'T98', 'T99', 'T101', 'T102', 'T103', 'T104', 'T105', 'T107', 'T108', 'T110', 'T113', 'T114', 'T117', 'T119', 'T120', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T131', 'T132', 'T133', 'T134', 'T136', 'T138', 'T139', 'T140', 'T141', 'T142', 'T143', 'T145', 'T146', 'T147', 'T148', 'T149', 'T150', 'T151', 'T152', 'T153', 'T154', 'T155', 'T156', 'T157', 'T159', 'T160', 'T161', 'T162', 'T164', 'T165', 'T166', 'T167', 'T168', 'T169', 'T170', 'T171', 'T172', 'T174', 'T175', 'T177', 'T179', 'T180', 'T181', 'T182', 'T183', 'T184', 'T185', 'T187', 'T188', 'T189', 'T190', 'T191', 'T192', 'T193', 'T194', 'T195', 'T196', 'T198', 'T200', 'T201', 'T202', 'T203', 'T204', 'T206', 'T207', 'T208', 'T210', 'T213', 'T215', 'T216', 'T219', 'T220', 'T222', 'T223', 'T224', 'T225', 'T226', 'T227', 'T228', 'T229', 'T230', 'T236', 'T238', 'T239', 'T241', 'T242', 'T243', 'T245', 'T246', 'T247', 'T249', 'T250', 'T255', 'T256', 'T259', 'T262']
def get_taxa_names(): return taxa_names