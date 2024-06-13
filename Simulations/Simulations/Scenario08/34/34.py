#!/usr/bin/env python
from numpy import *
data_1 = [
array([33.667112144950785, 33.096269682139166]),
array([31.822267195681693, 28.652606906920955, 28.299646862595207, 30.8114710963583, 27.74404482496262, 27.59728542451639, 27.44296628826886, 27.160659746629257, 27.794799506498467, 27.032142730725088]),
array([31.353274648213727, 30.757159154954515]),
array([31.53130013123448, 32.352535498808315]),
array([30.89251487916356]),
array([30.112477651917843, 31.07026142399085, 26.92568438680901, 26.897225737593097, 28.04691179971768, 28.031064856553925, 27.925748713987172, 27.95017274335756, 27.683226131938667]),
array([31.843333328910393, 29.477535734039343]),
array([30.934746532598815]),
array([27.278357611435073, 24.201689927643162, 25.591902334069463, 27.687927270989483, 25.19533248418996, 26.636223562811512, 26.789947451141646, 26.45264588941911, 23.24546544708552, 27.924035292826172, 25.81243020260223, 27.940362613739858, 17.693195113084037, 19.859770517056347, 16.110824694282165, 16.760552481581414, 17.50956289925867, 18.501572302810164, 17.824260601658764, 17.238323430966332, 20.165450072774153, 18.929474573427417, 16.707306240200452, 19.45740013356322, 17.309719788538057, 17.935925202244068, 13.746253652753776, 12.99683205517309]),
array([24.54483428384897, 27.011494166673305, 25.3050779283558, 23.319789178881575, 25.05221769916098, 24.619903244678405, 23.91014982915673]),
array([28.90715263459718, 27.145547825980085, 26.823364788574622, 26.917582195858568]),
array([29.130614311326056]),
array([29.228457665162164, 23.288474551283265, 24.465519915976177, 27.617916147297695, 27.41385615136871, 24.20811683051839, 27.931737195787086, 24.249092442149287, 23.536698728310526, 26.107913003558103, 27.81889190515143, 24.32220800024532, 24.051305116403412, 23.595744833264753, 25.456069102820273, 23.5286807256681, 22.38766144223036]),
array([28.903564667633745, 26.613468532854807, 27.41745969558847, 26.74576384155459, 25.752459064746546, 26.62584405097725, 26.29006394549444, 28.000072457367946]),
array([25.28286978664227, 25.447834072051663, 27.278222211073935, 27.506453570959035, 25.779530639959795, 24.1070477377707, 24.712939223415447, 25.552623395141957, 24.94832154863359, 22.127319151575456]),
array([28.13163093257587]),
array([28.48081563612756, 27.70106700079429, 26.773846415129416, 27.760312138624187, 27.20410327227743, 26.856104828421067, 26.719118912221457, 26.468923166648892, 27.502814071412665, 27.65512982406691]),
array([28.78525045324151, 28.085695014701585, 27.7895852233332]),
array([28.18659662983445, 24.741592350750807, 26.544836756484013, 25.286932916992473, 24.72971990087194, 25.440651643805708, 24.592228918846537, 27.373090370005798, 26.93352223574166]),
array([27.871594358129293]),
array([28.604293976735327, 27.737289780195827, 27.419969326136354, 26.84490719339492, 26.12979053630799, 25.4389143390393, 24.878796954706768, 25.386891853698202, 26.66207756565935, 27.53526887327419, 26.240815941203206]),
array([26.143597059211718, 27.62621596508559, 25.458184895747777, 26.487230565410936, 27.270074733765306, 27.21942965916673, 26.24249006197034, 23.7526058838311, 27.502330958505237, 25.95847851158809, 23.824208079828693, 23.835416255802272, 27.853665981046657, 24.222446558129725, 24.8237956913825, 26.774192173887037, 27.832494726807663, 22.095680117352583, 20.961439132268897, 21.6228167049412, 19.95838138857892, 19.664250502651218, 16.128049952457065, 20.07091695262753, 18.154891003214242, 19.86026523903, 16.045858769159246, 16.379810381874925, 17.91733863312613, 16.59392836739186, 18.572568773680075, 19.826105803777512, 18.984418945984917, 15.771152307910521, 15.62009900556165, 13.560025922075516, 13.106447612657682, 11.853415878955673, 12.831724602712008, 12.969883392243021, 13.078851077644812, 13.577230090267255, 12.019749372986341, 11.828485018388648, 11.833679618261309, 12.658692345724495, 13.450103785484348, 12.886019822956678, 12.732771428967375, 12.945630731914958]),
array([28.068512119498095, 27.858168033550637, 27.68734562233272, 27.568302199856436, 27.514591177078557]),
array([27.10776868448951, 26.955842450197615, 27.846512771258904]),
array([27.84089285596779]),
array([27.493929321590233, 23.76297024954168, 23.893287806610413, 25.417040569577733, 25.323540764851625, 26.737135310475526, 23.72784979087694, 26.263450008167712, 27.28528833512866, 25.582403094597666, 27.233813038049963, 26.40215184686753, 24.201015169005803, 27.149015370481372, 22.379511811563034, 20.750843547923314, 22.106127224340053, 16.13803133217, 16.14094861371662, 19.69208724299517, 18.31605457641306, 16.946692286759735, 16.890319207749435, 18.405649043724175, 18.597649604364324, 19.81473963882249, 19.42976706685239, 17.828073240923185, 19.96022963142393, 18.486720327085028, 16.81500153465568, 16.210087266290362, 16.952125505178884, 19.79748973094159, 15.95316949394087]),
array([23.49188852361262, 24.872866010665952, 24.13229068890811, 25.80576415792824, 23.38410898781387, 26.343313559000197, 25.02670714338102, 24.96749631687857, 26.62145480942275, 25.822797240171, 25.209363321772546, 23.76041997895949, 25.41698340115154, 26.803158834553933, 23.64162854261258, 26.580414738381766, 27.641234992557642, 25.19198798992319, 26.33156427699805]),
array([25.429144254304116, 26.399733573157235, 25.822229549813905, 24.952167942479914, 27.154791170156223, 24.84255533715583, 25.247609308339143, 26.562809628043684, 26.467321859247946]),
array([27.64287153808074, 27.06294224343878, 27.426093685795742, 26.640223314769557, 26.361589940470157, 26.495661963717883]),
array([26.281067932817475, 24.6693193011067, 25.626088840507055, 26.120020344066052, 27.261780447151654, 24.28028651762766, 25.79401379274566, 23.11791311154572, 27.160304639621497, 24.839794014015386, 25.364493226911588, 27.06471530145615, 26.21192913653367, 24.25425260627258, 25.5263583014003, 26.5680290751849, 24.83156590463492, 21.09478116693003, 22.438985204569164]),
array([25.75624283784587, 26.455122784097167, 26.889777088067657, 24.61169513790751, 27.077236041011798, 24.28974698307585, 26.743291910494158, 24.645928051527864, 26.38713428889334]),
array([25.89257496999363, 25.90011298206246, 23.93581404149749, 23.58712513321639, 24.209587543940966, 24.341676872966868, 24.152296572348412, 25.666020576833, 24.58539502348443, 26.996648537560226, 24.835393597130203, 24.678635822511414, 26.901331317293128, 27.02535588058702, 26.17401056211527, 26.962728395036013, 23.580061242555846, 23.626984257748582, 24.243765482926015, 23.905784992083248, 26.096520986739836, 25.22850475454924, 24.41421186593689, 24.137670657038136, 23.319062282888854, 24.447125867169447, 26.52576236474629, 26.17414061831693, 24.915478341251784, 26.743806389022854, 26.458139657337444, 22.027439540175926, 21.55336621382013, 22.35161704961976, 22.2042899101295, 19.053839301647145, 16.058145205458185, 20.110784229143093, 16.248907335642773, 17.266551005386443, 16.744126847200658, 15.979341628456172, 19.06891414640146, 17.840920775554302, 19.24103156228394, 16.163880264633033, 18.360357936633406, 18.628376736670923, 19.71140043899497, 19.235574255512173, 14.665863674711689, 14.479350466309775, 14.8776264480261, 13.68078784764045, 13.694559376027609]),
array([24.704407616800907, 25.85476942391788, 26.84274187478872, 24.617093586343746, 24.125180655284375, 22.65446378892331]),
array([25.69421847523483, 25.62960716175355, 26.305420404683524]),
array([25.51946993390521, 26.912390121658582, 25.592684693003687, 25.023728617181504, 25.12005873113637, 23.612081538003654, 26.810971203183563, 26.35222616750806, 26.316562736475237, 24.568650588658418, 26.507500208099472, 24.210363579382616, 17.05723972643788, 17.131031171104336, 16.48241110794175, 18.03809042058204, 20.283484563925096, 19.423123042964054]),
array([26.279346600525535, 25.739220744265396, 26.638618666205968, 26.263627662399415, 26.661209730310194, 23.901491207153192, 26.121853836512283]),
array([23.983709806989562, 25.017919834803884, 23.909637617825016, 25.651020470200727, 23.404295523735648, 16.572515268470323, 19.814898986181262, 19.97313110126567, 18.12071032763476, 17.621135722465944, 19.27518941422958, 16.516511426230572, 16.35613158694357, 18.889300746802448, 16.199169114019035]),
array([24.67131423125713, 23.349349559066997, 24.33420704696452, 25.524298973260034, 25.921379671581366, 23.427695208160696, 25.528431694320545, 25.878024579506626, 24.348375004991453, 24.52732642940937]),
array([25.709642108586333, 24.090767390062897, 24.838469799857368, 24.089377775696217, 25.639359362660343, 24.478834845892997]),
array([23.738526473948212, 24.360971357008236, 25.15626766115347, 24.019136865922786, 23.46841646713775, 23.119227624953766, 25.37447147761492, 25.592326164923247, 24.634854902670575]),
array([25.521862693418818, 25.723035694559297, 25.69439642633366, 25.876067546319923, 25.738407725728255, 25.674307655697653]),
array([23.895667043300666, 24.76959902832354, 22.458796031804532]),
array([23.77045211550683, 24.2267776826175, 23.73409740772654]),
array([23.274435125299142, 23.818107929658442, 23.27347274551362, 21.86567903203802, 22.47876935736158, 22.196327843326337]),
array([24.139654731713062, 24.23954096593795, 23.875912297528828, 24.24793484297767, 19.758512960961106, 20.30315977953393, 19.408468335439395, 16.542280060193537, 17.449674164746817, 17.92791655900649, 18.889406608158577, 16.88075781119243, 16.101541891857963, 18.705124260618557, 18.30077902192142, 14.020757377335109, 13.524260501045113, 12.913220728758384, 13.503892170174508, 13.17808395757947, 11.841876041830075, 13.646584606634168, 12.433166555061442, 12.682008147428808, 13.103019673940743]),
array([23.88389924399992, 23.72999719701739, 23.897773754633732]),
array([23.772584556028246, 23.345067664771893, 23.402812766791307, 23.60186576273392, 22.53156179651468, 22.174718880396426, 18.848461634128764, 19.30272096840147, 18.526678357157373, 18.79230717275247, 18.711785193349815, 19.3684794529651, 20.378165758318783, 19.441504229150528, 19.947730977561633]),
array([22.23008148931633]),
array([23.121524628138438, 23.04469432454177]),
array([22.794309192199695, 17.945822239337186, 19.1572874439597]),
array([22.71241352063824, 20.270546084710652, 19.78311376919985, 18.918095499612118, 17.130177591806838, 19.20706193347695, 19.368469336058123, 17.183512815269168, 20.231896370059108, 14.681105137445424, 13.25735022431655, 12.619006866214107, 12.760425091706447, 11.775820754428773, 12.990917208469588, 7.7797278633113915, 11.051864506525812, 8.422981153417728, 6.575519196297835, 4.5166213512808, 3.935515225609098, 3.123186826396131, 0.9923556875576971, 0.9111161314940582, 0.6643181079623668, 0.6074379236975851, 0.6935323624336779, 0.0]),
array([21.823054067912572]),
array([20.62036658133979, 20.236502227244344, 20.107300475334107, 20.377197601748055, 19.451487350927902]),
array([20.80388998998294, 21.498146622484676, 20.19205487593999, 18.33721568212023, 18.38248855833637, 18.45429196923436, 19.25386358468271, 19.995568709942447, 20.421388194696572, 18.589888511373594, 18.84686155373397, 20.352188794960124, 19.506886795551853, 17.818043532255142]),
array([19.172178617370015, 16.440876085567876, 14.03229175143646, 13.761095413058985, 13.19113572941001, 12.613689017483022, 12.795028481436898, 12.845009887849898, 9.977411123040334, 8.132553850885532, 11.176704882580452, 7.414249166182345, 5.399116121883498, 4.083933389175309, 2.996883479101471, 3.0362739761934314, 1.6233020513438083, 0.9244841912861821, 0.29332587286290995, 0.0]),
array([21.581660781182837]),
array([19.105970187681418, 18.237458735504298, 16.61607235657785, 16.737212790562594, 19.309665568748002, 17.571987692413284, 16.160262799383993, 16.752798132686127, 19.040615575144074]),
array([18.255919562278454]),
array([17.13381256460974, 19.46592102077915, 18.450580174933293]),
array([20.070489435830552, 19.882187538089315, 20.03563412156423, 19.549956239619405]),
array([19.382473532641633, 18.259296842751596, 17.568813787154685, 16.529199358785664, 18.09827690863189, 18.853511849466337, 18.15400161755905, 17.464314333478804, 17.994432883686184, 17.337024434437602, 16.68397311248694, 13.535615181418626, 11.81370322460128, 12.732917576271936, 13.34638048000997, 13.558834047139612, 13.038660395441442, 11.753180467412232, 12.657444530569395, 13.616141712446684, 12.880529874579643, 7.934963161419267, 9.037009951806159, 10.935947353274699, 11.607756241283827, 10.94297535936045, 9.130727570379332, 8.07997615484978, 10.311626369797398]),
array([18.650683205357772, 17.92774935156723, 17.026731484137834, 19.081516329598813, 17.330799223627285, 17.000809530093605, 18.31999137144138, 17.050885361347884, 18.079698327893794, 17.53181496163443, 16.896559646786912, 18.30676317828052, 16.45649141160573, 16.19160767897688, 14.418829206496216, 14.982039757227266, 14.04216046646755, 13.789567134585193, 13.684703125200123, 13.248316885393718, 13.258924899866347, 13.457098431391127, 13.79020653536336, 13.326012775807566, 13.264600999813783]),
array([18.76893052042417, 19.02099591844366]),
array([17.874740404597794, 18.426207907116776, 18.226337055215485, 18.401340334315975, 17.832673618090332, 18.211504832727197]),
array([18.39351018436575]),
array([16.320440241200327, 17.32662096676631, 16.268798404077312, 16.07719066878561, 16.135282431705757, 17.51194391931292, 17.948898566514433, 17.697546903731265, 14.97196383958219, 15.857980603026318, 14.365631192411765, 13.502834183695647, 13.416123574820233, 12.592622118094233, 13.147001273191758, 11.863578910089533, 11.851665178908267, 11.95510087477362, 11.759895844261175, 13.807581454417415, 13.242232363012334, 10.040416158372473, 8.289624123437065, 8.404754811537888, 8.257290632087194, 11.553303257172304, 11.189590473574105, 8.000118886491276, 8.394600649474816, 8.894511055005657, 10.385407136743783, 11.06028075253521, 8.757088938037914, 6.087610467240291, 7.228756798562626, 6.122504319786919, 4.650191107120207, 5.069512530418362, 4.793011013571458, 4.4583429480687755, 4.8690690704018715, 4.830647853539987, 3.7884020601495734, 2.895272627440557, 3.448428267089889, 2.835343269016186, 1.5562919554892425, 0.7842090287242172, 0.8131270113319298, 1.4406128032917875, 1.2720480660228506, 1.2509519779647775, 0.7175613590975488, 0.13405141456241243, 0.13600830170549794, 0.4563307064871911, 0.20789829337421828, 0.1160831951026228, 0.0]),
array([17.082020814314355, 17.752720401529064, 17.86176651845252, 17.80678858862714, 17.840129801947967, 16.997577883794893, 16.796937438949453, 16.89419336026382, 17.29728241782898, 18.05867819576485]),
array([17.597441909134837, 17.889662048853197]),
array([16.29050695831735, 16.251637491360697, 16.15299666958197, 16.3383524548042, 14.718293174419813, 14.839533842069786, 15.772233756876489, 13.916747403854304, 15.606526979092607, 13.414299458652314, 13.19578822157085, 11.644338659089804, 12.881682498491552, 13.176286220241687, 13.010342032860626, 11.876928149188872, 13.082584571562897, 12.749523073470801, 13.106075833774325, 13.46795785742542, 12.609003352658428, 13.275092867394552, 12.428103411956293, 13.573632929437514, 11.749662558447449, 10.957796021172902, 10.046042363331612, 11.192419743109058, 10.086242213463464]),
array([16.52182185773476, 16.80972177716805, 17.385227820142592, 17.28533783813826, 17.24891231359397, 17.146628673517128, 16.837385216394278, 16.041669895385407, 16.74292410988085, 16.919887524795374, 13.113319307613756, 12.864090351883794, 12.627254843819156, 13.20850140674278, 13.206582413850688, 12.554543357804953, 12.992631100527213, 12.861549438479859, 13.098543145438546]),
array([17.09199729269492, 16.897372170576947, 17.004297876339056, 17.19445632447738, 15.654963027430146, 14.613547565547815, 13.509661960614396, 12.65668713964726, 12.875884146674053, 12.427277022137185, 12.894208512674059, 13.143800634266672, 12.085921030047942, 12.513822677883901, 8.79853423374302, 10.512620182854718, 7.866851409002337, 10.950252116547237, 8.160461661520321, 10.78376381448699, 9.719842643572465, 7.009389175632186, 5.910733254387697, 3.520288657450617, 2.2240154155437066, 1.1189584049051011, 1.1275965420916971, 1.7627967442594517, 1.7737986004586284, 0.9132147942383346, 0.3962425168842308, 0.5936993324841332, 0.49708727953697945, 0.5998661556660705, 0.5235131135694646, 0.36538333361600794, 0.0]),
array([16.586131363498744, 16.90102213606464, 16.92167084083476, 16.27432471507184, 16.290854261902965]),
array([16.352482353730437, 16.212204796312598, 16.190042826538185, 16.123811351340336, 16.26177022943848, 16.884985989698933, 16.690717850912925, 16.508549951906755, 14.722561203289073, 14.614820228586462, 15.685025966063971, 14.375281582868617, 13.730345967246183, 13.73879386767001, 12.895755106750256, 13.31718477475478, 12.965586990071143, 13.379370245293172, 12.90509922319243, 12.62235397478271, 11.677089389729598, 12.143856996640505, 12.631672340838872, 13.388783261279144, 12.759763574742163, 12.026492404805175, 12.94748623240734, 12.835905303840931, 13.50605902058996, 12.481144831443682, 12.903333410611738, 13.735101207327219, 13.348391309968685, 12.93313995228342, 11.770004327727333, 12.959006656178675, 9.301160518366075, 10.217584290322382, 11.57923018908833, 10.816899689805492, 7.377495133228218, 10.849447930186228, 8.713738632171719, 10.365305513238075, 10.577315388992963, 7.989023256501629, 10.777033079723376, 10.849721855627353, 7.642057799853991, 9.602239734833834, 6.13937235874022, 7.09700197688909, 5.771207485424029, 5.015032446859352, 4.782740042298579, 4.215471977816378, 4.402919257493277, 3.6687293074035843, 3.2406458043380564, 3.5477183005353634, 3.3619582295903694, 2.9880155863182045, 0.8037366880613543, 1.6914847883928923, 1.6079507220809246, 1.1567000490558574, 1.2470968482524531, 1.2959419176604423, 1.1808250796651265, 1.3456962988809449, 0.8959717204831781, 0.8237978640161585, 0.6142812594598839, 0.38967744197293114, 0.3624337784896181, 0.4191470323690221, 0.0]),
array([12.803887694030342, 12.64488686377913, 13.11341157561358, 13.032016563925595, 11.829455228816823, 13.476284578066283, 12.55530312757733, 12.330267496860433, 11.748520649702714, 11.647870198615912, 11.83368944525382, 11.935360007320423, 10.619545597649338, 9.93585532899839, 10.645488217613995, 8.176746012980745, 7.580284469528807, 9.82326792618459, 7.543492239886258, 3.74990658580282, 4.997393254454083, 4.815917170559601, 4.8182301124859634, 3.09568545759782, 2.8437933965269977, 1.4110197588356164, 0.8825572615794379, 0.6720394119717193, 0.48718638572778494, 0.2968467063154512, 0.0]),
array([16.01188890839331, 15.648795661817644, 13.312083669403458, 13.491519095199001, 13.475369808308601, 12.281008914674238, 13.337757086377874, 12.58427581131447, 12.102084702290808, 13.419572444123974]),
array([14.570012286148895]),
array([14.77562787533998, 13.057114789978828, 12.714335849122994, 13.290687437030414, 11.874232989789165, 11.864145882341345, 11.841434696515371, 10.873742822289962, 8.785613505804688, 9.681568279382537, 9.564986691843925, 10.429627738398082, 7.746575456870595, 5.415797717318901, 0.8564626214878939, 0.21351995979005056, 0.2127240316818645, 0.0]),
array([12.94384011115543, 11.699327941463665, 13.275127809961202, 12.110408577711052, 13.092626674278678, 12.401577650391062, 7.81415744532174, 7.389469796686601, 7.892055126331892, 9.554145312265907, 10.415087977614345, 9.193173648917348, 10.229360651076409]),
array([11.744779337666042, 12.007124838574867, 13.002589488359076, 12.640734650418294, 13.59508627294957, 12.273289907656796]),
array([12.881858227278473, 12.911317354802176, 13.632518004909981, 12.661620491139516, 11.970673700297567, 13.215391212498316, 13.805292873033578, 8.040276089074867, 9.569118605029804, 7.965308958386398, 8.10499233555334]),
array([14.453230823210154, 14.975602335596562, 13.902798094902643, 12.894255687722113, 12.43204854847527, 12.877643845258413, 13.04867322989582, 12.954583992430623, 12.52638374559819, 13.011363783328846, 13.59086891993632, 12.26093884768909, 12.682221967823702, 12.185214070891492, 13.684635870284234, 12.692874704559824, 12.05249061339589, 12.64440992043235, 12.476421559112111, 13.47545100361571, 11.3085610752474, 8.845897533092698, 8.03203902513051, 10.537775063613733, 9.17274024871825, 8.558981381098368, 9.606450095119108, 10.664053277863996, 5.9812539673563645, 3.992365251131538, 4.628376340858894, 3.4674247876771322, 3.111688968052232, 2.3596675021520443, 0.7839975861132775, 1.7743197937279716, 1.4081036601012764, 1.5669637331738688, 1.0217597013417623, 1.1679284595423582, 1.7450503336450338, 1.5426417859119337, 0.7983205169061445, 0.5293822941958608, 0.4802654998854358, 0.2925318571605658, 0.0]),
array([14.426987453735467, 14.056286366049457, 12.040930646791805, 13.464508357228123, 11.99443106346333, 13.31702594003834, 11.739045802534758, 13.35222376835685, 13.195030621408518, 13.541295018649476, 12.518746126538941, 12.63553906002763, 12.94063988651186, 11.724146164194693, 12.780970534656646, 11.900752894803663, 13.739927800639792, 12.662923875844102, 12.789355281790435, 12.895284854672251, 12.685568379196036, 12.231380928303235, 8.362527123303916, 10.970964585866394, 7.625156244660907, 10.36950269165645, 7.643768647993433, 11.042918756890723, 9.081093850025475, 9.48747926222963, 9.998699741340978, 9.823187602752952, 9.792257270946969, 8.575420367837644, 10.280913076025652, 9.894675498329832, 9.796073118363006, 9.409483371374195, 11.353101864152531, 9.154221969839766, 5.73026262563213, 5.784484986566754, 5.720365970922648, 4.036611718116334, 3.8544269264265285, 4.29807010798881, 4.356907222207929, 3.9790532276156063, 2.9975521431523955, 3.4451216359616494, 2.0591438904128445, 2.248390661483654, 0.8109101480173366, 1.252498838501218, 1.0586470251276823, 1.1672670554440137, 0.8522726185352596, 1.3189861558075315, 1.4531013171035037, 0.652898035838792, 0.22715033861383604, 0.5904467955944979, 0.32298784699107835, 0.7460911895401684, 0.4013146419037963, 0.07179475014862347, 0.0]),
array([14.752454044799833, 14.19543977924975, 14.334803980606155, 12.145535577197244, 12.369136473927355, 11.898266102300905, 12.103294871657031, 12.658415421001457, 13.042040064385603, 11.93339993851345, 12.141321703817844, 11.864666886140242, 12.331913083930901, 12.916957492136412, 13.115810793227885, 12.094463306933168, 13.468961565842298, 13.009525135751298, 12.293261275038194, 13.145533530653033, 13.735472664682032]),
array([13.936495678760126, 11.969276647562515, 12.818511988695114, 12.632540589063021, 12.895158404951644, 13.0004812930044, 13.439149719388093, 13.08970268453714, 13.557026788854964, 12.970383254410338, 13.002119405715087, 13.683235013775748, 13.495464489137191, 13.057515476213204, 12.889723252167194, 12.433349292253288, 13.659184998557823, 13.6653704987427, 10.578758022340725, 11.576255347582768]),
array([12.510346476133579, 13.027439744311906, 13.599442126910507, 13.40004059944906, 13.369593372675629, 13.48874670371381, 12.801466020278417]),
array([11.974525126016548, 12.199955212587902, 11.968630549189763, 12.786200239977068, 13.671225845915421, 11.708833348235448, 8.552193319648543, 11.571389612028202, 7.011864725791427, 3.1076460898267344, 2.793361193267072, 3.1836609146548214, 3.257210232474903]),
array([12.769461607618748, 11.655030384935172, 11.753817807017592, 12.759679508046634, 9.28395696843069, 9.751215377921236]),
array([11.939835461235297, 11.7116503856361, 12.199844414736154, 12.43262637650454, 12.422183693956304, 12.221186072836103, 12.260526280173757, 12.293741433781843, 13.300287657917481, 7.523936424510635, 10.06574010472723, 8.537546867576308, 9.6720020145827, 9.330731527442396, 9.882117643211886, 10.668697915690155, 10.68560346552726, 8.910363448161867, 7.146779196096684]),
array([12.77067828894714, 12.091696818986511, 10.388846731657699, 8.402377328910376, 8.322300366166793, 8.490052344057283, 8.066043787287667]),
array([11.909419011680928, 12.682216096585563, 11.890118963479232, 12.359390068184775, 11.839647350235364, 12.563095942269705, 12.423889615113808, 11.751414522501612, 11.973042441535322, 12.149758152530817, 12.593307637704036, 12.27492734536133, 8.72759997020977, 10.088443611085125, 10.280279526734303, 10.7971727662029, 7.995987457203993, 11.464300788472809, 9.998887105252784, 10.40007681599437, 11.380992373338948, 7.5894807210970665, 10.049246681206036, 9.137489207596671, 11.621767381491313, 9.8540849347218, 10.203594092991416, 6.725534109027683]),
array([12.015799223718522, 12.733074327201294, 12.025773279997516, 12.818306994838432, 12.927348947767614, 7.907353472469431, 10.823229518899225, 11.140541352908947]),
array([12.14115697582667, 11.822444911992285, 12.669076311771883, 11.640155655502717, 12.322237443808467, 10.693458745902818, 11.130070998897205]),
array([12.112472355334068, 11.785851752380495, 11.656898677726904, 11.800022893269746]),
array([11.246630155951905]),
array([11.7653586034073, 10.748156843260753, 10.10135637567502, 10.594772950873491]),
array([12.018347802558596, 11.651076049764372, 11.949390876223344, 11.94604650804036, 11.923140764634818, 9.738774201300691, 9.38315261809383, 6.829268186525132, 6.328397127844483, 6.888541986760626, 5.410064767171952, 3.339435515885084, 1.2730702993714096, 0.8712663130026366, 1.334930572494463, 0.8067696664898617, 0.6412487343163, 0.6535572720888306, 0.28712092343191997, 0.0]),
array([12.047586126460537, 11.675822838258716, 10.127318886516532, 9.852236512518742, 10.181875382928608, 10.403561955552886]),
array([11.935303099031666, 11.689274603571572, 11.659655595587267, 11.686719880537275, 11.935734065485798, 11.90029932302186]),
array([11.664482868114744, 11.652846811142993, 11.784882995990959, 11.72191311421923, 9.62599661273889, 10.487047861411803]),
array([10.085390936156696]),
array([11.124332572904002, 10.53721111767049]),
array([10.873974968508726, 9.779817305328148, 10.859449666400213, 8.57920937502323, 10.527180894639569]),
array([8.075117711144573, 10.390668725438724, 10.551331534549336, 7.256900350009888, 6.177768976602595, 5.018199150848694]),
array([6.312613602172696, 0.0]),
array([10.598729348294494, 10.269990247817685]),
array([10.790565384385504, 10.877230603611585]),
array([4.881445439524538, 4.452057049617698]),
array([10.150613564950332, 8.143413517485337, 9.475889557291833, 8.752641195729957]),
array([8.756127515803351, 8.288693671540843, 10.166812319697705, 8.765183044259881]),
array([6.586563822544024, 5.305619292881807]),
array([8.63593404353801, 9.074968178210645, 9.423054233348696, 9.18840420489179, 9.058066284207532, 9.536096128735988, 8.447689521215333]),
array([9.41355632396706]),
array([8.675925332013804, 7.443647549031914, 7.306610842898021, 9.291447550431641, 7.464444705102114, 6.472217478733251, 5.077507219255778, 3.915401775406745, 3.931442506524344]),
array([8.20952696346372, 9.697058459747264, 9.725243694801476, 6.335797047498748, 7.228524352782351, 4.977266408823738, 4.16146250036051, 2.6341243451081615]),
array([9.096928491922657, 8.975567183879017, 8.985724812090554, 8.106015725359113, 7.745595678047993, 5.932983812200032, 6.080638063756181, 5.452707197478658, 7.196061229139511, 7.156369844135931, 4.7746771940024155, 5.227404898461077, 3.837047440799915, 3.768263098774276]),
array([9.15589889429342, 6.936095624737683]),
array([7.86899274286383, 6.522095318290495, 2.6618867428040573, 2.435536868471085, 0.7035649509642699, 0.23483580046174657, 0.0]),
array([8.18655591763611, 7.8404235137120155]),
array([8.765839434860771]),
array([8.146178918725607, 7.725021322446622, 8.452898874553814]),
array([7.694189672761875, 3.630710591857056, 3.380110909112174, 2.888456508793338]),
array([8.36789262163281]),
array([8.246636429304715, 7.585453025263123, 6.592954819603697, 5.443360827925722, 5.830056531829157, 4.580959975700374, 3.988124163424298, 3.398294246774662, 3.5110451553976154, 1.0001020959138889, 1.4109079733262415, 0.9600903219743968, 1.540589566112425, 0.07986375217174951, 0.0]),
array([8.070333606192678]),
array([7.9601319064578115, 7.998675829225374, 7.345739635192638, 5.449003591098222, 7.057041097508349, 6.765323934913484, 5.5952651425831545, 3.347413873595091]),
array([8.074433548412534, 5.914439537486295, 6.186613730131716]),
array([7.663317631134739, 6.93279282719236]),
array([4.358605518001255, 4.341805607761598, 3.965156561320804, 3.0455497871792545]),
array([7.6166113331494945]),
array([6.00858834684023, 3.1187289779543472]),
array([7.023075732805331, 5.797457484305525, 6.312293569009461]),
array([7.133177959043977, 6.042616553031594, 3.780616639630371, 3.608800190391185, 2.1585037262537536, 0.6722377137439266, 0.0]),
array([4.441461165832643, 4.2708832783979, 5.323664924143392, 4.528223133604367, 4.969479622014848, 4.3845300606273465, 4.890018406368116, 2.7150735761306213, 2.6686704424184757, 1.8951758477102392, 1.8490660656975226, 1.5171203384244287, 0.9185293862930891, 1.6810984827798638, 1.554964278916443, 1.0721035571559274, 0.23298565572875118, 0.6501519952704874, 0.6264497474817338, 0.7081376247122583, 0.4981103134630828, 0.6455453897420987, 0.643527399372339, 0.7666998997972839, 0.721497478032378, 0.0]),
array([6.823674189765374, 6.754704220157784, 6.059979570904383, 4.077228113027876, 1.9946621050944908, 1.4285260354796256, 1.6604329065278436, 0.6750152715285207, 0.4666602295626288, 0.5565117796373905, 0.613261953096488, 0.0]),
array([6.919700539480082, 7.106639014690112, 6.267444940086195, 6.75309833228671]),
array([6.110359505197844, 6.247214113390652, 3.7139650994628015, 5.133653463841445, 5.312647819369495, 4.105502275906139, 2.6808019988744367, 2.9375983391764016, 1.6767903566357523, 1.7048031397404777, 0.9885562655242278, 0.9646408351139876, 0.8639818014569333, 0.6823786880799458, 0.4495464751916973, 0.31899051760840896, 0.6103862391236691, 0.0]),
array([5.47359453926414, 5.983926336642098, 5.649214113735106, 4.898320932320259, 2.8441046688660743, 2.8670047365882962, 3.018725221436787, 3.1087994846706732]),
array([2.778749162111292, 1.3562324404834567, 1.4094520072141747, 0.5048490359158692, 0.0]),
array([4.479581733620782, 2.660929503718946, 2.670258802506651]),
array([5.474492932095146, 5.875036209250355, 5.958165514869302, 4.714446070312901, 5.0147419302785945, 4.448306989098688, 3.068934069010104, 2.7594675081419715, 1.658493265651147, 1.3049858663044378, 1.5409708958524782]),
array([5.1864726509984065, 5.272801344946017]),
array([2.4405536027360397, 1.2412908552096584, 0.23386708004478718, 0.5931190354618953, 0.5365107579156462, 0.0]),
array([3.655447296949723]),
array([5.216250933334731, 5.234559598977609]),
array([3.5584916702193707, 2.9366265171885213, 1.4533941675554738, 1.0997775519600128, 1.0609613251388543, 0.2807865072133846, 0.4822175960466618, 0.58390330799691, 0.0]),
array([3.4596585236895327, 2.7662076067766876, 2.203039223046706, 1.7156491895022308, 1.5718980707628616, 0.6345238791072485, 0.0942094418294379, 0.0]),
array([4.176765928196604, 4.390713444039259, 3.739877249795519, 3.4149853751096226, 3.5860325808035878]),
array([4.202449018647488]),
array([2.9991050043410294, 2.204879004132036]),
array([4.365046783644565, 4.593803808999879, 3.4390537014059057, 3.5852712164973024, 3.121741000046527, 3.4731199225276463, 2.9014996040476406]),
array([4.158829675984948, 4.317163549494783, 3.6944310408588565, 3.2642863192319473, 3.2255040386560143, 3.5202073894941797, 2.9337859428704487, 3.5013357164682324, 3.5492187411481555, 0.9469542625722597, 0.5404870606375538, 0.5545925682562047, 0.6970088087878387, 0.16231354140302767]),
array([3.743135376847012]),
array([2.652733631317999, 2.6311489441525513, 2.930064135632138]),
array([2.9770290081562294]),
array([2.4563497613937004]),
array([3.719074245661112]),
array([3.278281620407141]),
array([3.372574414616034, 0.7863529967014924, 1.6432545758807717, 1.0303656929010185, 1.1865685867353117, 0.5591741457716372, 0.3006169480851071, 0.4970652935231683, 0.24657204234760288, 0.6622421934714112, 0.6255766773714364, 0.7199577315449348, 0.4140169985322548, 0.25713354458318993, 0.2113186752493993, 0.1127575443697722, 0.0]),
array([3.4492397559185766, 3.3752295417411036, 2.581460917219405, 1.5994628208825066, 0.7391638520943167, 0.6718688093893679, 0.6874532110722062, 0.2563755178868713, 0.0]),
array([3.168776488927971, 1.3295432264289264, 1.2047948857114035]),
array([3.514783852219008, 2.723319508389313, 0.0]),
array([3.38278893748228, 0.3949107553294833, 0.26956268868706856, 0.0]),
array([2.6938315032028908, 0.4334218983756156, 0.20715378536541096, 0.0]),
array([3.211423824693915, 3.358801451283034]),
array([2.855290197259324, 0.944622284603115, 0.4620796052094983, 0.237250974112594, 0.0]),
array([3.3172123185180133, 2.1590690701352355]),
array([1.6718349131194254, 1.5352201913918198, 1.6300768844087843, 0.39728587329216086, 0.7039514760814816, 0.4341386175499299, 0.04868354767040414, 0.0]),
array([2.6062795108019876, 3.219512720123545, 3.1570968572232587, 2.797007125034035, 2.795816623171408, 1.6498514319233377, 1.1629171699515293, 1.272412726483001, 1.6489162612187411, 0.2577717976082753, 0.0]),
array([2.734428468845708]),
array([0.5911438899200191, 0.0]),
array([2.8919467168125124, 2.4624751850089495, 0.0]),
array([1.03667432470772, 1.690041560666194, 1.2650934487277694, 0.6012089050150181, 0.0]),
array([0.9166267583541129, 1.4875521470702298, 1.2921953190913293, 1.1398993029142601, 0.5292399104379069, 0.0]),
array([1.375923631460691, 1.3356656903911608, 0.5455998429120532, 0.0]),
array([0.2281323382454279, 0.0]),
array([1.910183975749709, 1.738648871924693, 1.4730447573540042, 1.6391929731874026, 1.2309594651418716, 1.3798246156687175]),
array([1.086664262186705, 0.368128866795435, 0.0]),
array([2.032969565030615, 0.84149491896608, 1.251334797563595, 0.8080084000090522, 1.0894919156625966, 0.29897663174301203, 0.5052845418224777, 0.4492687056816517, 0.0]),
array([2.0710616219709204]),
array([1.3966528562822562, 1.327817342029696]),
array([1.3178863764925663, 0.9397254891506884]),
array([1.487192329781717, 1.3005149362560402, 1.522401629587613, 0.293299306821998, 0.6584600992407941, 0.5787538484091233, 0.3631565423192945, 0.20340535735807597, 0.6923912978539548, 0.4815495465191803, 0.3425029474751064, 0.08615674926147666, 0.0]),
array([0.9913152610671129]),
array([1.3000930662799521, 1.191077595888105, 1.6302507262566825, 0.5960391292763643, 0.5905573480455657, 0.17934253622445318, 0.3534392225089992, 0.7697095317251136, 0.0]),
array([0.8296527400141755, 1.3546661077989066, 1.5696710007216819, 0.6579084420175978, 0.663945754536047, 0.4971885302854003, 0.5919929988922615, 0.09247135242029537, 0.0]),
array([0.9320935151391678, 1.2009368963542302, 1.4232059560788555, 1.0226452760190914, 1.203024620701021, 0.6905732857050535]),
array([0.9403273976124173, 1.3334644016126087, 0.5908028246460633, 0.2955559723259612, 0.0]),
array([1.2928255039073682, 0.18045923345257575, 0.5762247861662735, 0.7071936462885393, 0.4507801208994441, 0.0]),
array([1.2316042052689133, 0.0]),
array([1.1597700297059297]),
array([0.920491692281706, 0.6416979547748026, 0.6642702548022383, 0.24992808336188577, 0.1499019512095482, 0.0]),
array([1.0281139933782704, 0.30371717691546096, 0.0]),
array([0.8310456862251261, 0.2984095129195916, 0.6619961118781624, 0.3155364900037003, 0.4686153374237972, 0.0]),
array([0.5133255130699026, 0.23195818629024212, 0.18400377520126177, 0.2578067813031075, 0.0]),
array([0.4615132952193312, 0.14386980570470764, 0.09053575880258527, 0.0]),
array([0.6152814965584333]),
array([0.20463397677046002, 0.0]),
array([0.2054703226629844, 0.0010663736749905867, 0.0]),
array([0.13766347871179513, 0.0]),
array([0.02122656715051531, 0.0])
]
d = [data_1]
names = ["34"]
def get_data(i): return d[i]
def get_out_name(i): return names[i]
taxa_names = ['T0', 'T1', 'T2', 'T3', 'T5', 'T8', 'T9', 'T10', 'T11', 'T13', 'T14', 'T15', 'T16', 'T18', 'T21', 'T22', 'T24', 'T25', 'T27', 'T28', 'T29', 'T31', 'T32', 'T33', 'T34', 'T35', 'T36', 'T37', 'T38', 'T39', 'T40', 'T41', 'T42', 'T43', 'T44', 'T45', 'T47', 'T48', 'T49', 'T50', 'T51', 'T52', 'T53', 'T54', 'T55', 'T56', 'T57', 'T58', 'T60', 'T61', 'T63', 'T65', 'T66', 'T67', 'T68', 'T69', 'T71', 'T72', 'T73', 'T74', 'T75', 'T76', 'T77', 'T78', 'T79', 'T80', 'T81', 'T82', 'T83', 'T85', 'T87', 'T88', 'T89', 'T90', 'T93', 'T94', 'T95', 'T96', 'T97', 'T99', 'T100', 'T103', 'T104', 'T107', 'T109', 'T110', 'T111', 'T112', 'T113', 'T115', 'T116', 'T117', 'T119', 'T120', 'T121', 'T122', 'T123', 'T124', 'T125', 'T126', 'T127', 'T128', 'T129', 'T130', 'T133', 'T136', 'T137', 'T139', 'T140', 'T142', 'T144', 'T145', 'T146', 'T147', 'T149', 'T152', 'T153', 'T155', 'T156', 'T157', 'T158', 'T159', 'T161', 'T162', 'T163', 'T164', 'T165', 'T166', 'T167', 'T169', 'T170', 'T172', 'T173', 'T174', 'T175', 'T176', 'T177', 'T178', 'T180', 'T182', 'T183', 'T185', 'T187', 'T188', 'T189', 'T190', 'T192', 'T193', 'T194', 'T195', 'T196', 'T197', 'T201', 'T202', 'T203', 'T205', 'T206', 'T207', 'T208', 'T209', 'T210', 'T211', 'T212', 'T213', 'T214', 'T216', 'T217', 'T218', 'T219', 'T220', 'T221', 'T222', 'T223', 'T224', 'T225', 'T226', 'T227', 'T228', 'T229', 'T231', 'T232', 'T235', 'T236', 'T237', 'T238', 'T239', 'T240', 'T241', 'T243', 'T244', 'T245', 'T246', 'T247', 'T250', 'T251', 'T253', 'T256', 'T257', 'T258', 'T260']
def get_taxa_names(): return taxa_names