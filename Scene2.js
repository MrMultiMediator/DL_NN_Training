class Scene2 extends Phaser.Scene {
	constructor(){
		super("playGame")
	}

	create(){
		this.background = this.add.tileSprite(0, 0, config.width, config.height, 'background')
		this.background.setOrigin(0,0);

		//This.players is an array of "tuples" player objects. The 1st element in the tuple is
		//the player object containing all attributes. The 2nd element is the phaser physics object for that player
		this.players = [];

		//Initial try with number of hidden neurons. Setting to number lasers plus 1, for the player
		this.nHidden = gSettings.maxSpears + 1;

		//Number of input neurons
		this.nInput = 2 + gSettings.maxSpears*3;		

		for (var i=0; i < gSettings.nPlayers; i++){
			this.players.push([new player(this.nHidden), this.physics.add.image(10.0, config.height/2, "player")])
			this.players[i][1].setCollideWorldBounds(true);
		}

		this.text = this.add.text(20, 40, "", {
			font: "15px Arial",
			fill: 'White'
		})

		//Create a variable to listen for keyboard events
		this.cursorKeys = this.input.keyboard.createCursorKeys();

		//evol is an array of evolution data objects that will contain average survival time, average standard deviation per player,
		//and maximum average survival time for each generation
		this.evol = []

		this.time = 0;
		this.spears_ever = 0;
		this.se = 0
		this.generation = 0;
		this.show_weights = false; //Whether or not to print the weights, biases, and player info at the end of the generation

		//Array of "tuples" of spear objects which contain the data in the first element, and the
		//image in the second element.
		this.spears = []
	}

	update(){
		//Scroll the background
		this.background.tilePositionX += 7.0;

		//Gravitation aceleration
		this.gravity();

		this.time++;

		this.updatePlayerTime();

		//Potentially generate a spear to attack the player
		this.spearGen();

		//Control the player's motion
		this.movePlayer();

		this.spearDel();

		this.collisionDetect();

		this.checkRestart();

		if(this.cursorKeys.space.isDown){
			this.show_weights = true
		}
	}

	//Gravitational acceleration
	gravity(){
		for (var i = 0; i < this.players.length; i++){
			if (this.players[1][0].state == 'alive'){
				this.players[i][0].vy += 10;
			}
		}
	}

	//Control player motion
	movePlayer(){
		for (var i = 0; i < this.players.length; i++){
			if (this.players[i][0].state == 'alive'){

				//Construct the array of input neurons for propagation, using values
				//normalized such that they have approx. zero mean and std dev of 1.
				var input = [(this.players[i][1].y-300.)/300., (this.players[i][0].vy-500.)/500.]
				for (var j = 0; j < gSettings.maxSpears; j++){
					input.push((this.spears[j][1].y-300.)/300.);
					input.push((this.spears[j][1].x-400.)/400.);
					input.push((this.spears[j][0].vx+gSettings.spearAvgVel)/gSettings.spearVelStd);
				}

				//Propagate and refresh the output neuron
				this.players[i][0].NN.propagate(input,3*gSettings.maxSpears+2, this.nHidden);

				if (this.players[i][0].NN.output >= 0.5){
					this.players[i][0].vy -= 20.;
				}

				this.players[i][1].setVelocityY(this.players[i][0].vy);

				//If we hit the top or bottom, kill the player
				if (this.players[i][1].y >= 520.0 || this.players[i][1].y <= 85.0){
					this.die(i);
				}
			}
		}
	}

	//Update survival time for living players
	updatePlayerTime(){
		for (var i = 0; i < this.players.length; i++){
			if (this.players[i][0].state == 'alive'){
				this.players[i][0].stime[this.players[i][0].stime.length-1]++;
			}
		}
	}

	//Every P frames, there will be a probability, Q, that a spear will be generated with a
	//velocity sampled from a distribution. P,
	spearGen(){
		var P = gSettings.spearEvery
		var Q = gSettings.spearProb
		var nSpears = this.spears.length

		while (nSpears < gSettings.maxSpears){
			var spearY = 85+Math.random()*415;	//Make the y position of the spear a random location on the screen
			var spearX = 800;					//Set x position to far right of screen

			//Sample a Gaussian distribution for velocity
			var spearVx = -p5.prototype.randomGaussian(gSettings.spearAvgVel, gSettings.spearVelStd);

			//If the laser is too slow, or going in the opposite direction (which is possible with the Gaussian),
			//set velocity to the average
			if (spearVx > -5){
				spearVx = -gSettings.spearAvgVel;
			}

			//Add new spear
			this.spears.push([new spear(spearX, spearY, spearVx), this.physics.add.image(spearX, spearY, "spear")]);
			nSpears++;
			this.spears[nSpears-1][1].setVelocityX(this.spears[nSpears-1][0].vx);

			this.spears_ever++;

			//console.log('y: '+spearY);
			//console.log('vx: '+this.spears[nSpears-1][0].vx)
		}
	}

	//This function will remove spears, deleting the data in the array for that spear, 
	//and resizing the array. This will all be facilitated by the changing of the 
	//'todelete' boolean to true
	spearDel(){
		var nSpears = this.spears.length;
		for (var i = 0; i < nSpears; i++){
			//If the spear reaches the end of the screen, destroy it
			if (this.spears[i][1].x < -27.0 || this.alldead) {
				this.spears[i][0].todelete = true
			}
		}

		//Accomodate the fact that the array length is changing by decreasing the iterator variable,
		//i by 1 after removing a spear. Also break from the loop once the real end of the shortened
		//array is reached.
		for (var i = 0; i < nSpears; i++){
			if (i > this.spears.length - 1){
				break;
			}
			if (this.spears[i][0].todelete == true){
				this.spears[i][1].destroy(); //Destroy phaser physics object
				this.spears.splice(i,1); //Remove all data for this spear from the spears array
				i--;
			}
		}
	}

	//Kill players that get hit by lasers
	collisionDetect(){
		for (var j = 0; j < this.spears.length; j++){
			for (var i = 0; i < this.players.length; i++){
				if (this.players[i][0].state == 'alive'){
					if (this.spears[j][1].x < 43 && this.spears[j][1].y < this.players[i][1].y + 43 && this.spears[j][1].y > this.players[i][1].y - 43){
						this.die(i);
					}
				}
			}
		}
	}

	//Set player q's status to dead, and delete the coresponding phaser physics object
	die(q){
		this.players[q][0].state = 'dead';
		this.players[q][1].destroy();
		//console.log("player "+q+" survived for "+math.mean(this.players[q][0].stime)+" on average +/- "+math.std(this.players[q][0].stime));
	}

	checkRestart(){
		//Check if we need to create a new generation. If so, create one.
		this.alldead = false
		for (var i=0; i < gSettings.nPlayers; i++){
			if (this.players[i][0].state == 'alive'){
				break;
			}
			if (i==gSettings.nPlayers-1) this.alldead = true
		}

		//Resurrect the generation of players if all are dead and they haven't been exposed to 1500 lasers yet
		if (this.alldead && this.spears_ever < 500){
			console.log('\n'+this.spears_ever)
			console.log(this.spears_ever-this.se)
			console.log(this.players.length)
			console.log(this.time);
			this.se = this.spears_ever
			this.revive_all(1)
		}
		//New generation
		if (this.alldead && this.spears_ever >= 500){
			//console.log(this.show_weights)
			var std_scale = 0.5 //How many standard deviations to look ahead for player survival time cutoff

			//Construct array of average survival time and standard deviation for each player, and also get necessary
			//evolution data for this generation, storing it in an evol_data object that is stored in the 'evol' array.
			var temp = this.build_stime()
			this.evol.push(new evol_data(temp[0], temp[1], temp[2]))

			//Compute the average survival time over players
			//Compute the standard deviation in average survival time over players
			this.spears_ever = 0;
			this.generation++;

			//Create a new generation of players. The functions inside of this function take an array of players, as well
			//as an evol_data object, and use this information to evaluate fitness and create a new set of players
			var temp2 = this.new_gen(this.players, this.evol[this.evol.length-1], std_scale)
			this.players = []
			this.players = temp2
			this.revive_all(0);		//Bring back the phaser physics images, representing our new players.

			//If there are somehow too many players, remove the last one
			while (true){
				if (this.players.length <= gSettings.nPlayers) break;
				this.players = this.players.pop()
			}
			console.log(this.players)

			//Display all the code needed to be typed by the developer if she plans on killing the code and rerunning where she left off
			if (this.show_weights) this.disp_data(std_scale);

			console.log(this.evol)
			console.log(this.generation)
			this.show_weights = false;
		}
	}

	//Revive all players
	revive_all(stat){
		for (var i=0; i < gSettings.nPlayers; i++){
			this.players[i][0].state = 'alive';
			if (stat == 1) this.players[i][0].stime.push(0);
			if (stat == 0) this.players[i][0].stime = [0];
			this.players[i][0].vy = 0;
			this.players[i][1] = this.physics.add.image(10.0, config.height/2, "player")
		}
	}

	//This function computes average survival time and standard deviation for each player, and also computes the averages over players
	build_stime(){
		var average = 0.0;
		var std_dev_data = [];
		var stime = [];
		for (var i=0; i < gSettings.nPlayers; i++){
			this.players[i][0].st_ave = math.mean(this.players[i][0].stime);
			this.players[i][0].st_std = math.std(this.players[i][0].stime);
			stime.push(this.players[i][0].st_ave);
			std_dev_data.push(this.players[i][0].st_std);
		}
		//Compute the average standard deviation per player from its respective average survival time
		var std_dev = math.mean(std_dev_data);

		//Compute the average survival time per player as well as maximum average survival time
		average = math.mean(stime);
		var max_stime = math.max(stime);

		return [average, std_dev, max_stime]
	}

	//Create a new generation of players
	new_gen(players, evol, std_scale){
		var temp1 = this.clone(players, evol, std_scale);							//Make copies of select number of genetically fit players for new set of players
		var temp2 = this.mutate(players, evol, std_scale, temp1.length, 500);				//Make mutated versions of select number of genetically fit players for new set of players
		var temp3 = this.breed(players, evol, std_scale, temp1.length+temp2.length, 500);	//Breed most fit players to fill in the remaining number of most fit players

		return temp1.concat(temp2.concat(temp3))
	}

	//Clone some fit players. There is an analytic function that determines the probability that a previous player will be cloned exactly in next generation
	//Exit the function after one run through of every player.
	clone(players, evol, std_scale){
		console.log("Cloning procedure running.")
		var new_players = []

		//Generate array of player ID's (player array index number)
		var seed_players = this.select(players, evol, std_scale)

		//Construct list of new players by extracting the identified players from "seed_players" array
		for (var i = 0; i < seed_players.length; i++){
				//new_players.push(players[i])
				//new_players[i][1] = 0

				//Create a 'new' player, and then set the neural network to be the same as the player from which it is being mutated
				new_players.push([new player(this.nHidden), 0]);
				Object.assign({__proto__: new_players[new_players.length-1][0].NN.__proto__}, players[i][0].NN);
		}

		console.log("Cloning procedure complete.")
		return new_players;
	}

	//Mutate some fit players. There is an analytic function that determines the probability that a previous player will have mutated versions.
	//Exit the function after a certain percentage (set in the function) of the target number of players has been reached.
	//If there is only one fit player (unable to breed), exit the function after the target number of players has been reached.
	//If there are no fit players, take the best one and make extreme mutations to this player to create the entire set of players.
	//current = number of players currently in new generation
	mutate(players, evol, std_scale, current, target){
		//Once we reach the target number of total players in the new generation, exit the mutation procedure
		console.log("Mutation procedure running.")
		var new_players = []
		var seed_players = []

		//shortcount keeps track of how many times there have not been enough (at least 1) players selected for mutation.
		//If it gets past 5, it will select a random player for mutation, and reset the counter back to 0.
		var shortcount = 0

		//Generate a vector of player IDs called seed players: players selected for their fitness, and create a mutated version of each.
		//Repeat ad infinitum until the target number of players have been reached.
		while (true){
			seed_players = this.select(players, evol, std_scale);
			if (seed_players.length > 0) {
				var mmag = 1./seed_players.length;
			} else {
				var mmag = 1.
			}
			for (var i = 0; i < seed_players.length; i++){
				//Create a 'new' player, and then set the neural network to be the same as the player from which it is being mutated
				new_players.push([new player(this.nHidden), 0]);
				Object.assign({__proto__: new_players[new_players.length-1][0].NN.__proto__}, players[i][0].NN);

				//Mutation. The magnitude of the mutation is inversely proportional to the number of seed players. If there are only a few (or only one),
				//the mutation magnitude should be large to create a diverse set of children. If there are many seed players, the mutation magnitude need
				//not be that huge.
				new_players[new_players.length-1][0].NN.mutate(mmag, 3*gSettings.maxSpears+2, this.nHidden);
				shortcount = 0

				if (new_players.length + current >= target) break;
			}
			//Mutate a random player if the selection routine didn't find anybody suitable 5 times in a row
			if (seed_players.length == 0) shortcount++;
			if (shortcount >= 5){
				new_players.push([new player(this.nHidden), 0]);
				Object.assign({__proto__: new_players[new_players.length-1][0].NN.__proto__}, players[math.floor(math.random()*(gSettings.nPlayers-.01))][0].NN);
				new_players[new_players.length-1][0].NN.mutate(mmag, 3*gSettings.maxSpears+2, this.nHidden);
				shortcount = 0;
			}
			if (new_players.length + current >= target) break;
		}

		console.log('Mutation procedure complete.')
		return new_players;
	}

	//Breed some fit players. There is an analytic function that determines the probability that two previous fit players will be bred.
	//Exit the function after the target number of players has been reached.
	//current = number of players currently in new generation
	breed(players, evol, std_scale, current, target){
		console.log("Breeding procedure running.")
		var new_players = [];
		var seed_players = [];

		//shortcount keeps track of how many times there have not been enough (at least 2) players selected for breeding.
		//If it gets past 5, it assumes there aren't enough fit players and will default to the mutate procedure
		var shortcount = 0

		while (true){
			if (new_players.length + current >= target) break;
			seed_players = this.select(players, evol, std_scale);
			if (seed_players.length <= 1){
				shortcount++;
			} else {
				//Add breed functionality
			}
			//Not enough breedable specimens. Default to mutation procedure for the remaining new players and exit function.
			if (shortcount >= 5) {
				var temp = this.mutate(players, evol, std_scale, current+new_players.length, target);
				console.log("Breeding procedure complete.")
				return new_players.concat(temp)
			}
		}

		console.log("Breeding procedure complete.")
		return new_players;
	}

	select(players, evol, std_scale){
		var seed_players = []
		var a = 1.0;
		var shift = evol.st_ave + std_scale*evol.st_std //This will make x = 0 for the analytic function that determines the probability be set to the std cutoff distance (minimum lookout point)
		var factor1 = math.log(2.)/(evol.st_max-shift)
		var select_counter = 0
		var fit_counter = 0

		for (var i = 0; i < players.length; i++){
			var x = (players[i][0].st_ave-shift)*factor1*(evol.st_std/players[i][0].st_std)

			//Avoid division by zero for players with extremely small average standard deviation.
			if (players[i][0].st_std < 0.1) x = (players[i][0].st_ave-shift)*factor1*(evol.st_std/0.1)

			var stat = 1.25*math.random();
			//Add player i to seed_players if the selection function value is greater than the random number stat
			if (math.exp(a*x)-0.9 > stat){
				seed_players.push(i)
				console.log(i+'   '+math.exp(a*x)+'   '+stat);
				select_counter++;
			}
			//Keep track of the number of players that were qualified, so we can determine how many of the qualified players ultimately were selected
			if (math.exp(a*x)-0.9 > 0){
				fit_counter++;
			}
		}
		console.log(fit_counter+'   fit players')
		console.log(select_counter+'   fit players selected for current procedure.')

		return seed_players;
	}

	disp_data(std_scale){
		//Show everything about the best performing players at the end of the generation (if requested
		//with spacebar, so that the developer can kill the code and pickup where she left off later)
		var aveplusshift = this.evol[this.evol.length-1].st_ave + 0.5*2.355*this.evol[this.evol.length-1].st_std
		aveplusshift = this.evol[this.evol.length-1].st_ave + std_scale*this.evol[this.evol.length-1].st_std
		var fitplayers = 0

		console.log('WARNING!!! REMEMBER TO CHANGE gSettings.nPlayers DOWN BELOW TO THE NUMBER OF PLAYERS ABOVE FWHM')
		console.log('ALSO REMEMBER TO EXTRACT THE evol_data object information in order to create a generation of players from this data')
		console.log('this.players = [];')
		console.log('for (var i=0; i < gSettings.nPlayers; i++){')
		for (var i = 0; i < gSettings.nPlayers; i++){
			//Only include the players that survived longer than the FWHM upper cutoff of Gaussian survival time distribution
			if (this.players[i][0].st_ave >= aveplusshift){
				fitplayers++;
				console.log('	this.players.push([new player(this.nHidden), 0])')
				console.log('	this.players[i][0].st_ave = '+this.players[i][0].st_ave)
				console.log('	this.players[i][0].st_std = '+this.players[i][0].st_std)
				console.log('	this.players[i][0].NN.bias = ['+this.players[i][0].NN.bias+']')
				console.log('	this.players[i][0].NN.weights[0] = math.matrix('+this.players[i][0].NN.weights[0]+')')
				console.log('	this.players[i][0].NN.weights[1] = math.matrix('+this.players[i][0].NN.weights[1]+')')
			}
		}
		console.log('}')
		console.log(fitplayers+' fit players. Replace gSettings.nPlayers with this number in the loop definition above')
	}
}

class evol_data {
	//This class is for evolution data, providing info to the user regarding how the players are improving over time
	constructor(ave, std, maxx) {
		this.st_ave = ave;	//Average of average survival time over players in this generation
		this.st_std = std;	//Average standard deviation of survival time per player across all iterations
		this.st_max = maxx;	//Max average survival for current generation
	}
}

class player {
	//Player object has four attributes: 1. Whether it is alive or dead,
	//2.) Survival times (in frames) for each iteration stored in array,
	//3.) Average survival time  4.) Survival standard deviation
	//5.) Y velocity.			 6.) Associated neural network
	constructor(nHidden) {
		this.state = 'alive';
		this.stime = [0];
		this.st_ave = 0.
		this.st_std = 0.
		this.vy = 0;
		this.NN = new NeuralNet(3*gSettings.maxSpears+2, nHidden);
	}
}

class NeuralNet {
	activate(val){
		return 1./(1.+Math.exp(-val));
	}

	propagate(input, ni, nh){
		//This function will receive the input neuron values as an array, and
		//propagate them thru the neural net.
		var temp = 0.0
		var z = [];
		z.push(math.zeros(nh))

		//Fill in the activations for the input neurons with the values passed in as arguments
		for (var j=0; j < ni; j++){
			this.a[0] = math.subset(this.a[0], math.index(j), input[j]);
		}

		//Sum up weights and activations from input layer to get z values in hidden layer
		z[0] = math.multiply(this.weights[0], this.a[0])

		//Apply activation function to hidden layer neurons
		for (var k=0; k < nh; k++){
			temp = this.activate(math.subset(z[0], math.index(k))+this.bias[0]);
			this.a[1] = math.subset(this.a[1], math.index(k), temp)
		}

		z.push(math.zeros(1));
		z[1] = math.multiply(this.weights[1], this.a[1]);

		//Apply actiavtion function to output neuron
		this.output = this.activate(math.subset(z[1], math.index(0))+this.bias[1]);
	}

	mutate(mmag, ni, nh){
		//This function will mutate the neural net's weights and biases. The highest magnitude of the mutation is passed in as 'mmag' (mutation magnitude)
		//mmag ranges from 0 to 1.

		//Bias mutation
		var mutation = 0.0;

		mutation = this.bias[0]*mmag*(math.random()-0.5);
		this.bias[0] += mutation;
		mutation = this.bias[1]*mmag*(math.random()-0.5);
		this.bias[1] += mutation;

		//Weight mutation
		var orig = 0.0;
		for (var k=0; k < nh; k++){
			for (var j=0; j < ni; j++){
				orig = math.subset(this.weights[0], math.index(k,j));
				mutation = (math.subset(this.weights[0], math.index(k,j)))*mmag*(math.random()-0.5);
				this.weights[0] = math.subset(this.weights[0], math.index(k,j), orig+mutation);
			}
		}
		for (var j=0; j < nh; j++){
			orig = math.subset(this.weights[1], math.index(0,j));
			mutation = (math.subset(this.weights[1], math.index(0,j)))*mmag*(math.random()-0.5);
			this.weights[1] = math.subset(this.weights[1], math.index(0,j), orig+mutation)
		}
	}

	constructor(ni, nh) {
		this.weights = [] 	//Array of weight matrices
		this.bias = []		//Array of biases
		this.a = [];
		this.a.push(math.zeros(ni)) //Array of input neuron activations
		this.a.push(math.zeros(nh)) //Array of output neuron activations

		//Create the first weights matrix connecting hidden layer with input layer
		//Number of rows = nh (number of hidden neurons)
		//Number of columns = ni (number of input neurons)
		this.weights.push(math.zeros(nh,ni));
		this.bias.push(-1.);

		//Create the set of weights connected the input neurons to the hidden neurons
		for (var k=0; k < nh; k++){
			for (var j=0; j < ni; j++){
				this.weights[0] = math.subset(this.weights[0], math.index(k,j), Math.random());
			}
		}

		//Create the second weights matrix connecting the hidden neurons to the output neuron
		this.weights.push(math.zeros(1,nh))
		this.bias.push(-1.)

		//Create the set of weights connecting the hidden neurons to the output neurons
		for (var j=0; j < nh; j++){
			this.weights[1] = math.subset(this.weights[1], math.index(0,j), Math.random())
		}

		this.output = 0.0 //Temporary placeholder so that this quantity copies over with object.assign()
	}
}

class spear {
	constructor(x, y, vx) {
		this.x = x;
		this.y = y;
		this.vx = vx;
		this.todelete = false;
	}
}