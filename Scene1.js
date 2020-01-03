class Scene1 extends Phaser.Scene {
	constructor() {
		super("bootGame"); //'super' makes the class inhereit all the characteristics of the parent class
	}

	preload() {
		//"background" is the string that identifies the image
		this.load.image("background","assets/background.png")

		this.load.image("player","assets/player.png")

		this.load.image("spear","assets/spear.png")

		this.load.image("RESTART","assets/RESTART.png")
	}

	create() {
		this.scene.start("playGame")
	}
}