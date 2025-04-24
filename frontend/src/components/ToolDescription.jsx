import React from 'react';

export default function ToolDescription(props) {
    return (
      <div className={'container tool-description'}>
          <img className={'tool-reference-image'} src={props.imagePath} alt={'description'}/>
          <h3 className={'tool-title'}>{props.toolTitle}</h3>
          <p className={'tool-text-description'}>{props.textDescription}</p>
      </div>
    );
}