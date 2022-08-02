import pytorch_lightning as pl
import torch.nn as nn
import torch

class simple_dnn(pl.LightningModule):
    '''
    Standard 1 hidden layer dnn to test against (baseline)
    '''
    def __init__(self, input_size, output_size, layer_specs, optimizer,lr, criterion):
        super(simple_dnn, self).__init__()

        # linear regression
        if not layer_specs:
            self.model = torch.nn.Linear(input_size, output_size)

        else:
            layers = []

            # First layer with hidden activation
            layers.append(torch.nn.Linear(input_size,layer_specs[0]))
            layers.append(torch.nn.ReLU())

            # All the hidden layers given from the list 'layerSpecs'
            for i in range(len(layer_specs)-1):
                layers.append(torch.nn.Linear(layer_specs[i],layer_specs[i+1]))
                layers.append(torch.nn.ReLU())
            # Final layer with output size = 'output_size'
            layers.append(torch.nn.Linear(layer_specs[-1],output_size))

            # Package up model for easy forward function
            self.model = torch.nn.Sequential(*layers)

        if criterion == "l1":
            self.criterion = nn.L1Loss()
        elif criterion == 'mse':
            self.criterion = nn.MSELoss()
        else:
            raise Exception('Incorrect criterion. Choices are: "l1" or "mse"')

        self.lr = lr
        self.opt = optimizer.lower()
        self.output_size = output_size

    def forward(self, x):
        out = self.model(x)
        return out

    def eval_batch(self, batch):
        # Make predictions
        x, y, city, half_hour, set_point, standard_temp, charge = batch
        y_hat = self(x)

        # Evaluate predictions
        loss = self.criterion(y_hat, y)

        # Include the rest of the information so we can plot within "validation_epoch_end"
        return loss, city, y, y_hat, half_hour, set_point, standard_temp, charge

    def configure_optimizers(self):
        if self.opt == 'sgd':
            optimizer = torch.optim.SGD(self.parameters(), lr=self.lr)
        elif self.opt == 'adadelta':
            optimizer = torch.optim.Adadelta(self.parameters(), lr=self.lr)
        elif self.opt == 'adagrad':
            optimizer = torch.optim.Adagrad(self.parameters(), lr=self.lr)
        elif self.opt == 'adam':
            optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        elif self.opt == 'rmsprop':
            optimizer = torch.optim.RMSprop(self.parameters(), lr=self.lr)
        return optimizer

    def training_step(self, batch, batch_idx):
        # Evaluate batch
        loss,_,_,_,_,_,_,_ = self.eval_batch(batch)

        # Configure result
        self.log('train_loss', loss)
        return loss

    def validation_step(self, batch, batch_idx):
        # Evaluate batch
        val_loss, city, y, y_hat,  half_hour, set_point, standard_temp, charge = self.eval_batch(batch)

        # Configure result
        self.log('val_loss', val_loss)

        return val_loss, y, y_hat, city, half_hour, set_point, standard_temp, charge
